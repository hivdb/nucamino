NucAmino
========

NucAmino is a nucleotide to amino acid alignment program optimized for virus
gene sequences.


Compilation
-----------

NucAmino is a program written in [Go programming language][golang]. You need
to have Go installed to compile it. The installation of Go are varied in
different systems or even in same system. Therefore, we introduced
[Docker][docker] to simplify the building process. You don't need to have
a native Go installed if you have Docker.

The installation of Docker is very simple and it supports most modern
operating systems like Linux, MacOS and even Windows. Please visit
[its website][docker] to retrieve installation packages for your system.

Once you installed Docker, just type this single command and you'll have a
ready-for-use `nucamino` executable file under `./build` folder:

```bash
make build
```

Note: Windows users need [MinGW][mingw] to run `make` commands.

Download Binaries
-----------------

For convenience' sake, we provide pre-compiled executables for mainstream
systems. Those binary files can be found in release pages. The current
release page is [v0.1.1][latest].

[golang]: https://golang.org/
[docker]: https://www.docker.com/
[mingw]: http://www.mingw.org/
[latest]: https://github.com/hivdb/NucAmino/releases/tag/v0.1.1
