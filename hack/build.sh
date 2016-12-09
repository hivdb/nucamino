#! /bin/sh

set -e

cd /go/src/github.com/hivdb/nucamino
go get ./...
mkdir -p build
cd build
echo "`go build -v github.com/hivdb/nucamino/cmd/nucamino 2>&1` => ./build/nucamino"
