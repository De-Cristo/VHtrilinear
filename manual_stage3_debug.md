# Stage 3 Manual Debug Guide

This guide reproduces Stage 3 step by step without using `run_pipeline.sh`, so you can isolate whether a failure comes from the pipeline wrapper or from a specific MG5 / reweighting step.

The guide is split into two parts:
- **Outside the Singularity container**: commands you run on the host to enter the container and prepare files.
- **Inside `trilinear-boost.sif`**: the actual MG5 / reweighting commands.

Working directory:

```bash
/afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear
```

## A. Outside the Singularity container

These commands run on the host system first.

### 0. Enter the container

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear
source lhapdf_env.sh
apptainer exec --bind /afs --bind /cvmfs --bind /tmp --bind /eos \
  trilinear-boost.sif bash
```

After this command, you are **inside** the container shell.

## B. Inside `trilinear-boost.sif`

Everything below this point runs inside the container shell unless noted otherwise.

## 1. Prepare the environment

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear
```

If `python` is missing but `python3` exists, create a temporary shim before running MG5:

```bash
mkdir -p ~/bin
ln -sf "$(command -v python3)" ~/bin/python
export PATH="$HOME/bin:$PATH"
```

Quick checks:

```bash
command -v apptainer
command -v python
grep -n '^urlbase' lhapdf/bin/lhapdf
```

Expected:
- `apptainer` is on `PATH`
- `python` resolves to a usable interpreter
- `urlbase` ends with `/current/`

## 2. Generate `hz_MC`

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/MG5_aMC_v2_5_5
./bin/mg5_aMC < /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/cards/proc_mc_hz
```

Check:

```bash
test -d hz_MC && echo OK
```

Expected:
- `hz_MC/` exists
- MG5 finishes without errors

## 3. Patch `hz_MC` cards

Copy the run and parameter cards:

```bash
cp /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/cards/run_card.dat hz_MC/Cards/run_card.dat
cp /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/cards/param_card.dat hz_MC/Cards/param_card.dat
```

Set the quick test values:

```bash
sed -i "s|.*=.*nevents.*|  1000 = nevents ! Number of unweighted events requested|" hz_MC/Cards/run_card.dat
sed -i "s|.*=.*ebeam1.*|  6800 = ebeam1  ! beam 1 total energy in GeV|" hz_MC/Cards/run_card.dat
sed -i "s|.*=.*ebeam2.*|  6800 = ebeam2  ! beam 2 total energy in GeV|" hz_MC/Cards/run_card.dat
sed -i 's|False.*=.*store_rwgt_inf|True = store_rwgt_inf|' hz_MC/Cards/run_card.dat
sed -i "s|^lhapdf = .*|lhapdf = /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/lhapdf/bin/lhapdf-config #|" hz_MC/Cards/amcatnlo_configuration.txt
```

Check:

```bash
grep -n "nevents\|ebeam1\|ebeam2\|store_rwgt_inf" hz_MC/Cards/run_card.dat
grep -n '^lhapdf =' hz_MC/Cards/amcatnlo_configuration.txt
```

Expected:
- `nevents` is `1000`
- `ebeam1` and `ebeam2` are `6800`
- `store_rwgt_inf` is `True`
- `lhapdf` points to the local `lhapdf-config`

## 4. Run `gevirt.sh`

```bash
cp /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/trilinear-RW/gevirt.sh .
./gevirt.sh hz_MC
```

Check:

```bash
ls -l check_olp.inc proc_ml
```

Expected:
- `check_olp.inc` exists
- `proc_ml` exists

## 5. Generate `hz_ME`

Install the loop diagram filter:

```bash
cp /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/trilinear-RW/vvh-loop_diagram_generation.py madgraph/loop/loop_diagram_generation.py
```

Build the temporary process card:

```bash
cat > _proc_hz_me <<'EOF'
import model hhh-model-new
EOF
cat proc_ml >> _proc_hz_me
cat >> _proc_hz_me <<'EOF'
output hz_ME
quit
EOF
```

Run MG5:

```bash
./bin/mg5_aMC < _proc_hz_me
rm -f _proc_hz_me
```

Check:

```bash
test -d hz_ME && echo OK
```

Expected:
- `hz_ME/` exists
- MG5 finishes without errors

## 6. Build `check_OLP`

Go to the subprocess directory:

```bash
cd hz_ME/SubProcesses
```

Copy the reweighting code:

```bash
cp /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/trilinear-RW/makefile .
cp /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/trilinear-RW/check_OLP.f .
cp ../../check_olp.inc .
```

Detect the ME subprocess directory and copy the include files:

```bash
ME_SUBDIR=$(find . -maxdepth 1 -type d -name 'P0_*' | head -1)
cp "$ME_SUBDIR/pmass.inc" .
cp "$ME_SUBDIR/nsqso_born.inc" .
cp "$ME_SUBDIR/nsquaredSO.inc" .
```

Detect the MC subprocess directory and copy the remaining files:

```bash
MC_SUBDIR=$(find ../../hz_MC/SubProcesses -maxdepth 1 -type d -name 'P0_*' | head -1)
cp ../../hz_MC/SubProcesses/c_weight.inc .
cp "$MC_SUBDIR/nexternal.inc" .
```

Copy PDF libraries if needed:

```bash
cp -r ../../HHH-libs/libpdf.a ../../HHH-libs/Pdfdata ../../HHH-libs/libLHAPDF.a ../../HHH-libs/PDFsets ../lib/ 2>/dev/null || true
```

Build:

```bash
make OLP_static
make check_OLP
```

Check:

```bash
test -x check_OLP && echo OK
```

Expected:
- `check_OLP` exists and is executable

## 7. Pre-Step-7 Diagnostics

Run these checks before launching LO events. They help separate a launch syntax issue from a compile/environment issue.

Check the command files and PDF config:

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/MG5_aMC_v2_5_5
grep -n '^lhapdf =' input/mg5_configuration.txt
grep -n '^lhapdf =' hz_MC/Cards/amcatnlo_configuration.txt
```

Check the compiler and Python tools available inside the container:

```bash
command -v python
command -v g++
command -v gfortran
python --version
g++ --version | head -n 1
gfortran --version | head -n 1
```

Check that Boost headers exist:

```bash
test -f /usr/include/boost/shared_ptr.hpp && echo "BOOST headers present"
test -f /usr/include/boost/algorithm/string.hpp && echo "BOOST algorithm headers present"
```

Inspect the generated PDF wrapper that is failing during compile:

```bash
sed -n '1,80p' hz_MC/Source/PDF/pdf_lhapdf6.cc
```

Optional: if you want to see the exact compile command MG5 uses, rerun the launch once and stop at the failure. The debug log is usually written under:

```bash
ls hz_MC/run_*_debug.log
```

## 8. Generate LO events

Go back to the MG5 root:

```bash
cd ../../
```

Create the launch file:

```bash
cat > _launch_hz <<'EOF'
launch hz_MC
noshowerLO
done
quit
EOF
```

Run MG5:

```bash
./bin/mg5_aMC < _launch_hz
rm -f _launch_hz
```

Find the generated LHE file:

```bash
find hz_MC/Events -name 'events.lhe.gz' -o -name 'unweighted_events.lhe.gz'
```

Expected:
- One of the two `.lhe.gz` files is present

## 9. Manual fix if `pdf_lhapdf6.cc` is missing Boost includes

If Step 8 still fails during compilation in `hz_MC/Source`, patch the generated PDF wrapper before trying again.

Go to the generated PDF wrapper:

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/MG5_aMC_v2_5_5/hz_MC/Source/PDF
```

Back up the file:

```bash
cp pdf_lhapdf6.cc pdf_lhapdf6.cc.bak
```

Insert the Boost includes above `using namespace std;`:

```bash
sed -i '/^using namespace std;$/i #include <boost/shared_ptr.hpp>\n#include <boost/foreach.hpp>\n#include <boost/algorithm/string.hpp>' pdf_lhapdf6.cc
```

Verify the change:

```bash
sed -n '1,20p' pdf_lhapdf6.cc
```

Expected:
- `#include <boost/shared_ptr.hpp>`
- `#include <boost/foreach.hpp>`
- `#include <boost/algorithm/string.hpp>`

Then rerun the launch step:

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/MG5_aMC_v2_5_5
./bin/mg5_aMC < _launch_hz
```

If it still fails, inspect the debug log:

```bash
ls hz_MC/run_*_debug.log
tail -n 80 hz_MC/run_01_LO_tag_1_debug.log
```

## 10. Reweight events

Unzip the LHE file into the ME subprocess directory:

```bash
gunzip -c "$(find hz_MC/Events -name 'events.lhe.gz' -o -name 'unweighted_events.lhe.gz' | head -1)" > hz_ME/SubProcesses/events.lhe
```

Run the reweighter:

```bash
cd hz_ME/SubProcesses
./check_OLP
```

Check:

```bash
ls -l events_rwgt.lhe
```

Expected:
- `events_rwgt.lhe` exists

## 11. Fix `pdlabel` in `check_OLP.f` if you see `Unimplemented distribution= lhapdf`

If `./check_OLP` stops with:

```text
Unimplemented distribution= lhapdf
```

patch the generated reweighting driver to use the named PDF set that this code expects.

Edit `hz_ME/SubProcesses/check_OLP.f` and change:

```fortran
        pdlabel='lhapdf'
```

to:

```fortran
        pdlabel='nn23nlo'
```

You can do it with:

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/MG5_aMC_v2_5_5/hz_ME/SubProcesses
sed -i "s/pdlabel='lhapdf'/pdlabel='nn23nlo'/" check_OLP.f
make check_OLP
```

Then rerun:

```bash
./check_OLP
```

This should allow the driver to recognize the PDF distribution.

## 12. Copy outputs

```bash
cd /afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear
mkdir -p output
cp MG5_aMC_v2_5_5/hz_ME/SubProcesses/events.lhe output/
cp MG5_aMC_v2_5_5/hz_ME/SubProcesses/events_rwgt.lhe output/
ls -lh output/events.lhe output/events_rwgt.lhe
```

Expected:
- Both files exist in `output/`

## How to interpret failures

- If `hz_MC` generation fails, the issue is in the process generation or environment setup.
- If `gevirt.sh` fails, the issue is in subprocess extraction or MG5 config.
- If `hz_ME` generation fails, the issue is likely the loop-generation filter or model import.
- If `check_OLP` fails to build, the issue is in copied include files, PDF libs, or the makefile.
- If event generation fails, the issue is in MG5 launch syntax, `run_card.dat`, or LHAPDF/PDF access.
- If reweighting fails, the issue is in `check_OLP` inputs or the generated LHE file.
- The pipeline scripts now apply the Boost-header patch and the `pdlabel='nn23nlo'` fix automatically, so this manual section is primarily for debugging older runs or verifying the exact fallback steps.

## Host-side notes

If the container command fails before you get the shell prompt, the problem is outside the pipeline itself:

- `apptainer` socket or permission errors mean the container cannot start in the current environment.
- If MG5 complains about `/usr/bin/env python`, then the host/container environment is missing a `python` executable, even if `python3` exists.
