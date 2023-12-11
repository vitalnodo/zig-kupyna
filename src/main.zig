const std = @import("std");
const testing = std.testing;
const assert = std.debug.assert;
const tables = @import("tables.zig");

pub fn Kupyna(comptime n: u10) type {
    comptime assert(n >= 8 and n <= 512 and n % 8 == 0);
    return struct {
        const Self = @This();
        pub const digest_length = n / 8;
        pub const block_length = if (n <= 256) 512 / 8 else 1024 / 8;
        pub const rounds = if (n <= 256) 10 else 14;
        const columns = if (n <= 256) 8 else 16;
        const rows = 8;
        pub const Options = struct {};
        const STATE = [columns][rows]u8;
        state: STATE,
        buf: [block_length]u8 = undefined,
        buf_len: u8 = 0,
        total_len: u64 = 0,

        pub fn init(options: Options) Self {
            _ = options;
            if (block_length == 512 / 8) {
                var state: STATE = undefined;
                std.mem.writeInt(
                    u512,
                    std.mem.asBytes(&state),
                    1 << 510,
                    .big,
                );
                return Self{
                    .state = state,
                };
            } else {
                var state: STATE = undefined;
                std.mem.writeInt(
                    u1024,
                    std.mem.asBytes(&state),
                    1 << 1023,
                    .big,
                );
                return Self{
                    .state = state,
                };
            }
        }

        pub fn update(d: *Self, b: []const u8) void {
            var off: usize = 0;
            if (d.buf_len != 0 and d.buf_len + b.len >= block_length) {
                off += block_length - d.buf_len;
                @memcpy(d.buf[d.buf_len..][0..off], b[0..off]);
                d.round(&d.buf);
                d.buf_len = 0;
            }
            while (off + block_length <= b.len) : (off += block_length) {
                d.round(b[off..][0..block_length]);
                d.total_len += block_length;
            }
            const b_slice = b[off..];
            @memcpy(d.buf[d.buf_len..][0..b_slice.len], b_slice);
            d.buf_len += @as(u8, @intCast(b[off..].len));
        }

        pub fn hash(b: []const u8, out: *[digest_length]u8, options: Options) void {
            var d = Self.init(options);
            d.update(b);
            d.final(out);
        }

        fn round(d: *Self, b: *const [block_length]u8) void {
            var input: STATE = undefined;
            @memcpy(std.mem.asBytes(&input), b);
            var after_T_xor: STATE = undefined;
            for (0..columns) |j| {
                for (0..rows) |i| {
                    after_T_xor[j][i] = d.state[j][i] ^ input[j][i];
                }
            }
            T_xor(&after_T_xor);
            var after_T_add: STATE = input;
            T_add(&after_T_add);
            for (0..columns) |j| {
                for (0..rows) |i| {
                    d.state[j][i] ^= after_T_xor[j][i] ^ after_T_add[j][i];
                }
            }
        }

        fn count_bits(arr: []const u8) u14 {
            var i: usize = arr.len - 1;
            while (i < arr.len) : (i -%= 1) {
                const octet = arr[i];
                if (octet != 0) {
                    var bit: u8 = 7;
                    while (bit < 8) : (bit -%= 1) {
                        if ((octet & std.math.shl(u8, 1, bit)) != 0) {
                            return @intCast(i * 8 + bit + 1);
                        }
                    }
                }
            }
            return 0;
        }

        pub fn final(d: *Self, out: *[digest_length]u8) void {
            // const in = "404142434445464748494a4b4c4d4e4f5051000000000000000000000000000000000000000000000000000000000000000000008f0200000000000000000000";
            // _ = std.fmt.hexToBytes(&d.buf, in) catch unreachable;
            // @memset(d.buf[d.buf_len..], 0);
            // d.buf[d.buf_len - 1] = 0x51;
            const last_block_bits_count = count_bits(&d.buf);
            const total_bits = d.total_len * 8 + last_block_bits_count;
            const extra_bits = total_bits % 8;
            if (block_length - d.buf_len < 12) {
                d.round(&d.buf);
                @memset(d.buf[0..], 0);
            }
            if (extra_bits > 0) {
                const mask = ~std.math.shr(u8, 0xFF, extra_bits);
                const pad_bit = std.math.shl(u8, 1, (7 - extra_bits));
                d.buf[d.buf_len - 1] &= mask;
                d.buf[d.buf_len - 1] |= pad_bit;
            } else {
                d.buf[d.buf_len] = 0x80;
                d.buf_len += 1;
            }
            std.mem.writeInt(
                u96,
                d.buf[block_length - 12 ..],
                total_bits,
                .little,
            );
            d.round(&d.buf);
            var state = d.state;
            T_xor(&state);
            for (0..columns) |j| {
                for (0..rows) |i| {
                    d.state[j][i] ^= state[j][i];
                }
            }
            const state_as_bytes = std.mem.asBytes(&d.state);
            @memcpy(out, state_as_bytes[state_as_bytes.len - digest_length ..]);
        }

        fn T_xor(state: *STATE) void {
            for (0..rounds) |r| {
                kappa_v(r, state);
                pi(state);
                tau(state);
                psi(state);
            }
        }

        fn T_add(state: *STATE) void {
            for (0..rounds) |r| {
                eta_v(r, state);
                pi(state);
                tau(state);
                psi(state);
            }
        }

        fn kappa_v(v: usize, state: *STATE) void {
            for (0..columns) |j| {
                state[j][0] ^= @intCast((j << 4) ^ v);
            }
        }

        fn eta_v(v: usize, state: *STATE) void {
            var s: *[columns]u64 = @ptrCast(@alignCast(state));
            for (0..columns) |j| {
                var tmp: u64 = ((((columns - j - 1) * 0x10) ^ v) << (7 * 8));
                tmp ^= 0x00f0f0f0f0f0f0f3;
                tmp +%= s[j];
                s[j] = tmp;
            }
            _ = &s;
        }

        fn pi(state: *STATE) void {
            for (0..columns) |j| {
                for (0..rows) |i| {
                    state[j][i] = tables.sboxes[i % 4][state[j][i]];
                }
            }
        }

        fn tau(state: *STATE) void {
            for (0..rows - 1) |i| {
                _tau(state, i, i);
            }
            const shift = if (block_length == 512 / 8) 7 else 11;
            _tau(state, rows - 1, shift);
        }

        fn _tau(state: *STATE, row: usize, shift: usize) void {
            var tmp: [columns]u8 = undefined;
            var i: usize = columns - shift;
            while (i < columns) {
                tmp[i] = state[i][row];
                i += 1;
            }
            i = columns - shift;
            while (i >= 1) {
                state[(i - 1 + shift) % columns][row] = state[i - 1][row];
                i -= 1;
            }
            i = 0;
            while (i < shift) {
                state[i][row] = tmp[i + columns - shift];
                i += 1;
            }
        }

        fn psi(state: *STATE) void {
            var tmp: [rows]u8 = undefined;
            for (0..columns) |j| {
                @memset(&tmp, 0);
                var i: usize = rows;
                while (i >= 1) {
                    var product: u8 = 0;
                    var b: usize = rows;
                    while (b >= 1) {
                        product ^= tables.mult_table[state[j][b - 1]][tables.mds_matrix[i - 1][b - 1]];
                        b -= 1;
                    }
                    i -= 1;
                    tmp[i] = product;
                }
                for (0..rows) |_i| {
                    state[j][_i] = tmp[_i];
                }
            }
        }
    };
}

const h = "000102030405060708090A0B0C0D0E0F101112131415161718191A1B1C1D1E1F202122232425262728292A2B2C2D2E2F303132333435363738393A3B3C3D3E3F404142434445464748494A4B4C4D4E4F505152535455565758595A5B5C5D5E5F606162636465666768696A6B6C6D6E6F707172737475767778797A7B7C7D7E7F808182838485868788898A8B8C8D8E8F909192939495969798999A9B9C9D9E9FA0A1A2A3A4A5A6A7A8A9AAABACADAEAFB0B1B2B3B4B5B6B7B8B9BABBBCBDBEBFC0C1C2C3C4C5C6C7C8C9CACBCCCDCECFD0D1D2D3D4D5D6D7D8D9DADBDCDDDEDFE0E1E2E3E4E5E6E7E8E9EAEBECEDEEEFF0F1F2F3F4F5F6F7F8F9FAFBFCFDFEFF";
test "Kupyna256" {
    const test_input = &[_][]const u8{
        h[0 .. 512 / 8 * 2],
        h[0 .. 1024 / 8 * 2],
        h[0 .. 2048 / 8 * 2],
        "FF",
        h[0 .. 760 / 8 * 2],
        "",
        h[0 .. (510 / 8) * 2] ++ "3C",
        h[0 .. (655 / 8) * 2] ++ "50",
    };
    const test_output = &[_][]const u8{
        "08F4EE6F1BE6903B324C4E27990CB24EF69DD58DBE84813EE0A52F6631239875",
        "0A9474E645A7D25E255E9E89FFF42EC7EB31349007059284F0B182E452BDA882",
        "D305A32B963D149DC765F68594505D4077024F836C1BF03806E1624CE176C08F",
        "EA7677CA4526555680441C117982EA14059EA6D0D7124D6ECDB3DEEC49E890F4",
        "1075C8B0CB910F116BDA5FA1F19C29CF8ECC75CAFF7208BA2994B68FC56E8D16",
        "CD5101D1CCDF0D1D1F4ADA56E888CD724CA1A0838A3521E7131D4FB78D0F5EB6",
        "875C0023DAA0C077809FDD6A9672B49E03903BFF98EBE48740AE998C7BE3851E",
        "4237D7DE1A00C4CC8037EDE9C54BA60D1C705CD1495DE19E5245BF3509DB59CE",
    };
    for (0..test_input.len) |i| {
        var input: [256]u8 = undefined;
        const input_slice = try std.fmt.hexToBytes(&input, test_input[i]);
        var output: [256]u8 = undefined;
        const output_slice = try std.fmt.hexToBytes(&output, test_output[i]);
        const Kupyna256 = Kupyna(256);
        var k = Kupyna256.init(.{});
        k.update(input_slice);
        var result: [Kupyna256.digest_length]u8 = undefined;
        k.final(&result);
        try testing.expectEqualSlices(u8, output_slice, &result);
    }
}

test "Kupyna48" {
    const Kupyna48 = Kupyna(48);
    var input: [64]u8 = undefined;
    _ = try std.fmt.hexToBytes(&input, h[0 .. 512 / 8 * 2]);
    var out: [6]u8 = undefined;
    Kupyna48.hash(&input, &out, .{});
    const expected = [6]u8{ 0x2F, 0x66, 0x31, 0x23, 0x98, 0x75 };
    try testing.expectEqualSlices(u8, &expected, &out);
}
