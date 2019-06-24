# Containers

[Singularity](https://sylabs.io/singularity/) is a container solution, specifically tailored for HPC clusters.

Because of the many dependencies of `libomogrid` and `gnina`, a Singularity container allows to quickly install and use the software on different platforms (from local machines, to HPC clusters).

## Installation

### Singularity 2.4


## Build a Container

A Singularity container can be easily built from a single definition (`.def`) file:
```
singularity build CONTAINER.img CONTAINER.def
```

The build process requires `sudo` privilegies and might take a while. Once the container is built, it can be copied or moved on different platforms.

## Use

The `--nv` option to Singularity commands enables [Nvidia](https://www.nvidia.com/en-us/) support.

### Interactive Session

An interactive session is useful for compilation and installation purposes and can be lunched as wollows:
```
singularity shell CONTAINER
```

### Run Scripts

Scripts can be easily run inside a Singularity container as follows:
```
singularity exec CONTAINER ./SCRIPT
```