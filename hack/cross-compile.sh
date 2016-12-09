#! /bin/sh

set -e

cd /go/src/github.com/hivdb/nucamino
go get ./...
mkdir -p build
cd build
for GOOS in darwin linux windows; do
    for GOARCH in 386 amd64; do
        echo "`go build -v -o nucamino-${GOOS}-${GOARCH} github.com/hivdb/nucamino/cmd/nucamino 2>&1` => ./build/nucamino-${GOOS}-${GOARCH}"
    done
done
