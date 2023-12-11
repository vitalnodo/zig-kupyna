# kupyna hash

A **failed** attempt to implement the Kupyna hash function because I found that it 
uses a complex enough padding, which other implementations don't implement. The
reference implementation uses an explicit message bit length. Kupyna256 tests 
for N=760 and N=510 give incorrect results. Maybe something also doesn't work.

```bash
zig test src/main.zig
```

## links
https://eprint.iacr.org/2015/885.pdf

https://github.com/Roman-Oliynykov/Kupyna-reference