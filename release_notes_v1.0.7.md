# ROHan v1.0.7

A build fix. No heterozygosity estimate changes, and no behavior changes at run time: if
v1.0.6 built for you, this release does nothing you need.

## Fixed: the link failed with "undefined reference to libdeflate_crc32"

On a machine with `libdeflate-dev` installed, the build compiled everything and then failed at
the final link with a page of errors like:

```
/usr/bin/ld: ../lib/htslib/libhts.a(bgzf.o): in function `hts_crc32':
bgzf.c:575: undefined reference to `libdeflate_crc32'
/usr/bin/ld: ../lib/htslib/libhts.a(cram_io.o): in function `libdeflate_deflate':
cram/cram_io.c:872: undefined reference to `libdeflate_alloc_compressor'
```

htslib's `configure` looks for libdeflate and uses it whenever it finds it, so `libhts.a` ends
up referencing `libdeflate_crc32` and friends. ROHan's makefile, however, listed its link
libraries by hand, and that list had no `-ldeflate` in it. The build therefore broke on
exactly those machines that happened to have the development package installed, and nowhere
else — which is why it went unnoticed for so long and why it could not be reproduced on a
machine that did not have it.

Simply adding `-ldeflate` would have moved the problem rather than fixed it: machines without
libdeflate would then fail with `cannot find -ldeflate`. The flag is needed sometimes and
wrong at other times, and which it is depends on what htslib's configure found.

htslib already writes out the exact set of libraries its static build needs, into
`htslib_static.mk`, and samtools consumes that file for precisely this reason. ROHan now does
the same, so the link line follows whatever htslib was actually configured with — libdeflate
or not — instead of a hand-maintained guess.

**If you worked around this** by adding `-ldeflate` to `src/Makefile` yourself, you can drop
the change; it is no longer needed and no longer harmful either way.

**bam2prof needs the same treatment** and is not covered by this release, as it lives in a
separate repository. If your build of it hit the same wall, the `-ldeflate` workaround there
is still the answer for now.

## Also

* The installation instructions said `cd ROHan` after a `git clone` that creates `rohan`.
  Fixed.

## Downloading

The attached `rohan-v1.0.7-linux-x86_64.tar.gz` is a statically linked Linux x86-64 build.
Unpack it and run `bin/rohan` from inside the unpacked directory:

```bash
tar xzf rohan-v1.0.7-linux-x86_64.tar.gz
cd rohan-v1.0.7-linux-x86_64
./bin/rohan --help
```

ROHan locates its deamination, error, base-frequency and pre-computed coverage files relative
to the binary, so keep the directory intact rather than copying `bin/rohan` out on its own.

To build from source instead:

```bash
git pull
make
make test
```

**Thanks** to the reporter for a write-up that included the full linker output and the gcc
version — that is what made the cause obvious.
