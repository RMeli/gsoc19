# Tools


## C++ Programming

### GDB

[GDB](https://www.gnu.org/software/gdb/)

[GDB-CUDA]()

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

### Valgrind

[Valgrind](http://valgrind.org/) is a system for debugging and profiling Linux programs which can automatically detect many memory management and threading bugs.


## Scripting

### GNU Parallel

[GNU Parallel](https://www.gnu.org/software/parallel/) is a shell tool for executing jobs in parallel. A job can be a single command or a small script that has to be run for each of the lines in the input.

The following bash function
```
BASH_FUNCTION(){
    echo $1
}

# Export function for GNU parallel
export -f BASH_FUNCTION 
```
can be easily run in parallel on `NUM_CPUS` CPUs for a list of different `INPUTS`
```
parallel -j NUM_CPUS BASH_FUNCTION ::: INPUTS
```