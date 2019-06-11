# Compilation and Installation

## OpenBabel

### Download

OpenBabel needs a few patches in order to work properly. Use the following fork: `https://github.com/RMeli/openbabel/tree/fix/icode`.

```
git clone https://github.com/RMeli/openbabel
cd openbabel
git checkout fix/icode
```

### Compilation

```
mkdir build && cd build
cmake .. \
    -DWITH_JSON:BOOLEAN=FALSE \
    -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DPYTHON_EXECUTABLE=/usr/bin/python3 \
    -DCMAKE_INSTALL_PREFIX=$HOME/software/openbabel
make -j
ctest
```

### Installation

```
mkdir -p $HOME/software/openbabel
make install
```

The following variables have to be set in order to compile `smina` (and other software) with this custom installation of OpenBabel:
```
export LIBRARY_PATH=$HOME/software/openbabel/lib
export LD_LIBRARY_PATH=$HOME/software/openbabel/lib
```

`LIBRARY_PATH` is used before compilation to include directories containing static and/or shared libraries that have to be linked to the program. `LD_LIBRARY_PATH` is used by the program (successfully compiled and linked) to dinamically link shared libraries.

## Smina

## CMake