#/bin/sh

# please run these interatively
exit 1

# go to the Tc_c_tools directory and build the c code
pushd ../..
make
popd


# install rdkit
pip install rdkit
pip install numpy

# a command line inferace package
pip install fire
