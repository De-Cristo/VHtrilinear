# Rebuild Local Container and Rerun Stage 3

This guide rebuilds `trilinear-boost.sif` with Boost headers installed, then reruns Stage 3.

**When to use:** `hz_MC/Source/PDF/pdf_lhapdf6.cc` fails to compile because `boost/shared_ptr.hpp` is missing.

> **Note:** This does NOT push to Docker Hub. You build everything locally and output a `trilinear-boost.sif` that drops in as a replacement. You do not need Charlotte's Docker Hub account.

> **Note:** The [`Dockerfile`](file:///afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/trilinear-env/Dockerfile) already contains `boost-devel boost-static` on the `yum install` line — no edits needed unless you want to add more packages.

---

## 1. Check the Dockerfile (already updated)

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/trilinear-env
grep "yum -y install" Dockerfile
```

Expected output:
```
RUN yum -y install gcc-gfortran boost-devel boost-static
```

If the output is missing `boost-devel`, edit `Dockerfile` to add it. Otherwise proceed.

---

## 2. Build the Docker image locally with podman

On lxplus, Docker is not available but `podman` is. Use a writable tmp directory to avoid storage driver conflicts:

```bash
mkdir -p /tmp/$USER/podman-root /tmp/$USER/podman-runroot /tmp/$USER/podman-tmp
export TMPDIR=/tmp/$USER/podman-tmp
export PODMAN_TMPDIR=/tmp/$USER/podman-tmp

podman --root /tmp/$USER/podman-root \
       --runroot /tmp/$USER/podman-runroot \
       --storage-driver=vfs \
       build --isolation chroot \
       -t trilinear-local:latest \
       -f Dockerfile .
```

Expected: build completes with `Successfully tagged trilinear-local:latest`.

**If you get a storage database mismatch** (`overlay does not match vfs`), reset the tmp dirs and retry:
```bash
rm -rf /tmp/$USER/podman-root /tmp/$USER/podman-runroot /tmp/$USER/podman-tmp
# then repeat the mkdir + export + build commands above
```

---

## 3. Export the image and build the SIF

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling

# Export docker image to a tar archive
podman --root /tmp/$USER/podman-root \
       --runroot /tmp/$USER/podman-runroot \
       --storage-driver=vfs \
       save -o /tmp/$USER/trilinear-local.tar trilinear-local:latest

# Convert to SIF (drop-in replacement — local filename)
apptainer build trilinear-boost.sif docker-archive:///tmp/$USER/trilinear-local.tar
```

Expected: `trilinear-boost.sif` appears in the selfcoupling directory.

---

## 4. Verify Boost is inside the new SIF

```bash
apptainer exec trilinear-boost.sif \
    bash -lc 'test -f /usr/include/boost/shared_ptr.hpp && echo BOOST_OK || echo BOOST_MISSING'
```

Expected: `BOOST_OK`

Also confirm gfortran is present:
```bash
apptainer exec trilinear-boost.sif bash -lc 'gfortran --version'
```

---

## 5. Re-run Stage 3

With the new SIF in place, re-run `run_pipeline.sh` (it will pick up `trilinear-boost.sif` automatically):

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear
./run_pipeline.sh --nevents 1000    # quick smoke test
```

Or follow the manual debug guide for step-by-step control:
- [`manual_stage3_debug.md`](file:///afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/manual_stage3_debug.md)

---

## 6. Success criteria

- `pdf_lhapdf6.cc` compiles without Boost errors
- `./bin/mg5_aMC < _launch_hz` completes and produces `events.lhe.gz` in `hz_MC/Events/run_01_LO/`
- `./check_OLP` produces `events_rwgt.lhe`
- `output/events.lhe` and `output/events_rwgt.lhe` appear

---

## 7. Cleanup (optional)

Once the SIF is built and verified, you can free up podman storage:

```bash
rm -rf /tmp/$USER/podman-root /tmp/$USER/podman-runroot /tmp/$USER/podman-tmp
```
