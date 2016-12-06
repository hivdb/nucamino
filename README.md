NucAmino
========

NucAmino is a nucleotide to amino acid alignment program optimized for virus
gene sequences.


Compilation
-----------

NucAmino is a program written in [Go programming language][golang]. You need
to have Go installed to compile it.

For users of Unix-like systems (Linux, MacOS, FreeBSD, etc), we recommend
installing the latest version of Go with [GVM][gvm]. Please access [its Github
page][gvm] to find instruction for your system.

For Windows users, Go official website provide [MSI installation files][gowin].
With [MinGW][mingw] installed, you can also install Go with GVM.

When you have the environment prepared, use this command to retrieve latest
code:

```bash
go get github.com/hivdb/nucamino/...
```

Then you can `cd` into the nucamino folder:

```
cd ${GOPATH}/src/github.com/hivdb/nucamino
```


[golang]: https://golang.org/
[gvm]: https://github.com/moovweb/gvm
[gowin]: https://golang.org/doc/install#windows
[mingw]: http://www.mingw.org/
