#! /bin/sh

set -e
BUILDFOLDER=${1:-"build"}
if [ -z $GOPATH ]; then
    GOPATH=${HOME}/go  # Enable default GOPATH (as per https://rakyll.org/default-gopath/)
fi

cd ${GOPATH}/src/github.com/hivdb/nucamino
go get ./...
mkdir -p $BUILDFOLDER
cd $BUILDFOLDER
echo "`go build -gcflags="-l=4" -compiler="gc" -v github.com/hivdb/nucamino 2>&1` => ./$BUILDFOLDER/nucamino"
