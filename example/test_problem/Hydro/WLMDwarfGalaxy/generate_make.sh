# This script should run in the same directory as configure.py

PYTHON=python3

${PYTHON} configure.py --mpi=true --hdf5=true --fftw=FFTW3 --gpu=true --model=HYDRO --mhd=true \
                       --particle=true --gravity=true --flu_scheme=MHM_RP --flux=HLLD --passive=2 \
                       --par_attribute_flt=2 --dual=ENPY --star_formation=true --feedback=true --grackle=true \
                       --bitwise_reproducibility=true --libyt=true --nlevel=12 "$@"
