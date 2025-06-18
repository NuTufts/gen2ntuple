#!/bin/bash

GEN2NTUPLE_DIR=$PWD
export GEN2NTUPLE_BINDIR=${GEN2NTUPLE_DIR}/build/installed/bin
export GEN2NTUPLE_LIBDIR=${GEN2NTUPLE_DIR}/build/installed/lib
echo "Gen2ntuple directory: ${GEN2NTUPLE_DIR}"

export PATH=${GEN2NTUPLE_BINDIR}:${PATH}
export LD_LIBRARY_PATH=${GEN2NTUPLE_LIBDIR}:${LD_LIBRARY_PATH}

