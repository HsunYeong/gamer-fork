# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --model=HYDRO --flu_scheme=MHM_RP --flux=HLLD --mhd=true --slope=PLM \
                       --cr_streaming=true --eos=GAMMA --hdf5=true --openmp=true --double=true "$@"
