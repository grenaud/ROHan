# Cutting a release

Two separate things have to be statically linked and neither happens by default:

* `bin/rohan`, which is tracked in git;
* the **release tarball**, which people download instead of cloning.

Both matter because ROHan is run on clusters where libgsl is absent from the compute nodes,
and a dynamically linked binary fails there with `libgsl.so.NN: cannot open shared object
file`. The ordinary `make` produces a dynamic binary, so the static link has to be done
deliberately before tagging. That is the step that is easy to forget and hard to notice
afterwards.

The worked example below is v1.0.5; substitute the version you are releasing.

## 1. Start from a clean, pushed master

```bash
git checkout master
git pull
git status --porcelain | grep -v '^??'      # must print nothing
```

## 2. Build everything

```bash
make
```

The first build on a fresh clone clones and builds htslib, samtools, libgab and libharu into
`lib/`, so it needs network access and takes a while. This leaves `bin/rohan` **dynamically**
linked; the next step replaces it.

## 3. Build the binary that ships

```bash
make -C src/ static
file bin/rohan
```

`file` must say `statically linked`. Expect one linker warning:

```
warning: Using 'getpwuid' in statically linked applications requires at runtime
the shared libraries from the glibc version used for linking
```

It is benign — it concerns `getFullPath()` resolving `~`.

## 4. Test the binary you are about to ship

```bash
make -C testData/ test-clean
make -C testData/ test
```

Expect `all 20 checks passed`. A few minutes, no network needed.

Use `make -C testData/ test`, not the top-level `make test`: the top-level target depends on
the phony `bin/rohan` and runs `make -C src/` first. That does not relink in practice, because
`bin/rohan` is newer than every object file by this point, but it will relink — quietly
replacing the static binary with a dynamic one — if any source was touched after step 3.

## 5. Commit the binary

```bash
file bin/rohan                              # statically linked, one last time
git add bin/rohan
git commit -m "static build for v1.0.5"
git push
```

Check `file` again rather than trusting step 3. Anything that ran `make` in between undid it.

## 6. Build the release tarball

The binary cannot ship on its own. It resolves its default data files **relative to its own
location** (`getFullPath(cwdProg+"../...")`, `rohan.cpp` ~5174), so a bare binary fails at
startup or, worse, part way through a run:

```
ERROR: Cannot find file .../bin/../preComputated/coverageprior/correctionCov_8.bin
       containing pre-computed coverage weights
```

Four directories have to travel with it, laid out so that `../` from `bin/` finds them:

```bash
V=v1.0.5
D=rohan-$V-linux-x86_64
rm -rf /tmp/$D && mkdir -p /tmp/$D/bin
cp bin/rohan /tmp/$D/bin/
strip /tmp/$D/bin/rohan
cp -r deaminationProfile illuminaProf DNAprof preComputated /tmp/$D/
tar czf /tmp/$D.tar.gz -C /tmp $D
sha256sum /tmp/$D.tar.gz
```

`strip` halves the binary (8.7 MB to 4.4 MB) and the whole tarball comes to about 2 MB, since
the four data directories are only 260 KB together. Do not strip `bin/rohan` in the repository,
only the copy going into the tarball.

## 7. Smoke test the tarball somewhere else

Do not skip this. It is the only check that the packaged layout actually resolves, and it is
how the missing `preComputated/` above was caught:

```bash
rm -rf /tmp/unpack && mkdir /tmp/unpack && tar xzf /tmp/$D.tar.gz -C /tmp/unpack
cd /tmp/unpack
./$D/bin/rohan --size 500000 --chains 2000 --seed 42 \
   --auto  <repo>/testData/simulated.autosomes \
   -o /tmp/unpack/t \
   <repo>/testData/simulated.fa <repo>/testData/simulated.bam
```

It must run to completion and write `/tmp/unpack/t.summary.txt`. Run `<repo>/bin/rohan --help`
from the unpacked copy too and confirm the `--deam5p`, `--err` and `--base` defaults point
inside the unpacked directory rather than at your build tree.

## 8. Tag, push the tag, publish

```bash
git tag -a v1.0.5 -m "ROHan v1.0.5"
git push origin v1.0.5
```

A tag is not pushed by a plain `git push`; it needs the second command.

Then create the release and attach the tarball:

```bash
gh release create v1.0.5 /tmp/rohan-v1.0.5-linux-x86_64.tar.gz \
   --title "ROHan v1.0.5" \
   --notes-file release_notes_v1.0.5.md
```

Paste the sha256 from step 6 into the notes so people can verify the download. The same can be
done through the GitHub web interface if `gh` is not available.

To confirm afterwards that a given fix is really in a release:

```bash
git tag --contains <commit>
```

Empty output means the commit is on master but in no tagged release — worth checking before
telling anyone to "upgrade" rather than "build from master".

## Checklist

- [ ] master clean and pushed
- [ ] `make -C src/ static` run
- [ ] `file bin/rohan` says statically linked
- [ ] `make -C testData/ test` passes 20/20
- [ ] `file bin/rohan` still says statically linked
- [ ] binary committed and pushed
- [ ] tarball built with all four data directories, stripped binary inside
- [ ] tarball unpacked elsewhere and run to completion
- [ ] sha256 recorded
- [ ] tag created **and** pushed
- [ ] release published with the tarball attached
