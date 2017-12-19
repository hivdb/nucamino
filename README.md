NucAmino
========

[![Travis CI](https://api.travis-ci.org/hivdb/nucamino.svg?branch=master)](https://travis-ci.org/hivdb/nucamino)
[![codecov](https://codecov.io/gh/hivdb/nucamino/branch/master/graph/badge.svg)](https://codecov.io/gh/hivdb/nucamino)

NucAmino is a nucleotide to amino acid alignment program optimized for virus
gene sequences.

The paper: https://www.ncbi.nlm.nih.gov/pubmed/28249562


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
systems. These binary files can be found in release pages. The current
release page is [v0.1.3][latest].

Usage
-----

Once compiled or installed, you can use NucAmino to analyze virus sequences
by using a command line program "nucamino". For now, we only support HIV
pol sequences but more may be added in the future. The instruction of the
command line tool can be retrieved by the following command:

Linux and MacOS:

```shell
./nucamino --help
./nucamino hiv1b --help
```

Windows:

```windows
nucamino --help
nucamino hiv1b --help
```

[golang]: https://golang.org/
[docker]: https://www.docker.com/
[mingw]: http://www.mingw.org/
[latest]: https://github.com/hivdb/NucAmino/releases/tag/v0.1.3
