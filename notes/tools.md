# Tools

## C++ Programming

### GDB

[GDB](https://www.gnu.org/software/gdb/)

[GDB-CUDA](https://developer.nvidia.com/cuda-gdb)

#### Run Program

Hook program to GDB:

```bash
gdb PROGRAM
```

To run the program, use `run`:

```bash
(gdb) run
```

To hook a program to GDB, including its command-line arguments:

```bash
gdb --args PROGRAM ARGS
```

The program can be run by simply typing `run` as before.

#### Set Breakpoint

Set a breakpoint in `FILE` at line number `LINE`:

```bash
(gdb) break FILE:LINE
```

#### Remove Breakpoint

Information about all breakpoints can be accessed via

```bash
info break
```

Single breakpoints can be removed as follows:

```bash
del BREAKPOINT_NUMBER
```

#### Explore Variables

Print variable value:

```bash
(gdb) print VARIABLE
```

Print variable type:

```bash
(gdb) ptype VARIABLE
```

#### Continue Execution

In order to resume execution after hitting a breakpoint, use

```bash
(gdb) continue
```

Sometimes it's useful to skip loops or other constructs. In order to achieve this, the execution can be restarted until a given line

```bash
(gdb) until LINE_NUMBER
```

#### Debug a Segmentation Fault

When a `Segmentation Fault` occurs within GDB, the backtrace can be printed with

```bash
(gdb) backtrace
```

or

```bash
(gdb) where
```

### Valgrind

[Valgrind](http://valgrind.org/) is a system for debugging and profiling Linux programs which can automatically detect many memory management and threading bugs.

## Scripting

### GNU Parallel

[GNU Parallel](https://www.gnu.org/software/parallel/) is a shell tool for executing jobs in parallel. A job can be a single command or a small script that has to be run for each of the lines in the input.

The following bash function

```bash
BASH_FUNCTION(){
    echo $1
}

# Export function for GNU parallel
export -f BASH_FUNCTION
```

can be easily run in parallel on `NUM_CPUS` CPUs for a list of different `INPUTS`

```bash
parallel -j NUM_CPUS BASH_FUNCTION ::: INPUTS
```
