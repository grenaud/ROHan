# ROHan v1.0.8

A build fix for bam2prof, which ROHan's top-level `make` downloads and builds. The `rohan` source is unchanged since v1.0.7, so heterozygosity and ROH estimates do not change.

## Fixed: bam2prof failed to link with "undefined reference to libdeflate_*"

The v1.0.7 notes listed this as still open. On a machine with `libdeflate-dev` installed, ROHan built `bin/rohan`, but the build then failed linking `bam2prof/bam2prof` with the same `libdeflate_crc32` / `libdeflate_alloc_compressor` errors that v1.0.7 fixed for rohan itself. bam2prof's own Makefile lists its link libraries by hand, with no `-ldeflate`, but here it is linked against the htslib that ROHan configured.

The top-level Makefile now gives bam2prof the libraries recorded in htslib's `htslib_static.mk`, as `src/Makefile` already does for rohan. The link line therefore matches how htslib was actually configured, with or without libdeflate. If you added `-ldeflate` to `bam2prof/Makefile` by hand, you can drop that change.

## Also

* bam2prof is fetched with `git clone --depth 1` (retried up to three times) instead of `wget` of the GitHub zipball, which can return 502/504 errors for minutes at a time. `unzip` is no longer needed.
* The README has a "How to cite" section with BibTeX and badges.

## Downloading

The attached `rohan-v1.0.8-linux-x86_64.tar.gz` is a statically linked Linux x86-64 build, against htslib 1.25 and samtools 1.25. Unpack it and run `bin/rohan` from inside the unpacked directory:

```bash
tar xzf rohan-v1.0.8-linux-x86_64.tar.gz
cd rohan-v1.0.8-linux-x86_64
./bin/rohan --help
```

ROHan finds its deamination, error, base-frequency and pre-computed coverage files relative to the binary, so keep the directory intact rather than copying `bin/rohan` out on its own.

sha256: `0d94b9b1105769ff15caa2295a6fda221ca25e7654c88c4f25926b267ecb43b7  rohan-v1.0.8-linux-x86_64.tar.gz`

To build from source instead:

```bash
git pull
make
make test
```
