# Compilation and Installation


## OpenBabel

[Open Babel](http://openbabel.org/wiki/Main_Page) is a chemical toolbox designed to speak the many languages of chemical data. It's an open, collaborative project allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state materials, biochemistry, or related areas.

### Requirements

* [CMake](https://cmake.org/)

Optional:

* [Eigen3](http://eigen.tuxfamily.org/)
* [Python](https://www.python.org/)

```
git git build-essential \
    cmake \
    libeigen3-dev \
    python3-dev
```

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

[smina](https://sourceforge.net/projects/smina/) is a fork of [Autodock Vina](http://vina.scripps.edu/) that focuses on improving scoring and minimization.

### Requirements

* [Boost](https://www.boost.org/)
* [Eigen3](http://eigen.tuxfamily.org/)
* [Open Babel](http://openbabel.org/wiki/Main_Page)

```
apt install git build-essential libboost-all-dev libeigen3-dev
```

Install OpenBabel as described above.

### Download

```
git clone https://git.code.sf.net/p/smina/code smina
```

### Compilation

```
cd smina/build/linux/release
make -j
```

## CMake
