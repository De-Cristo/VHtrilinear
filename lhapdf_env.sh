#!/bin/bash
# Source this to set up LHAPDF paths, then optionally enter the container.
export PYTHONPATH="/afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/lhapdf/lib64/python2.7/site-packages/:${PYTHONPATH:-}"
export LD_LIBRARY_PATH="/afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/lhapdf/lib/:${LD_LIBRARY_PATH:-}"
export PATH="/afs/cern.ch/user/l/lichengz/works_dir/private/selfcoupling/VHtrilinear/lhapdf/bin:$PATH"
