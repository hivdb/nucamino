#! /bin/sh

set -e
BUILDFOLDER=${1:-"build"} 

cd /go/src/github.com/hivdb/nucamino
go get ./...
mkdir -p $BUILDFOLDER
cd $BUILDFOLDER
echo "`go build -gcflags "-l=4" -v github.com/hivdb/nucamino/cmd/nucamino 2>&1` => ./$BUILDFOLDER/nucamino"
