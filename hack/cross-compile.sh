#! /bin/sh

set -e

cd /go/src/github.com/hivdb/nucamino
go get ./...
mkdir -p build
cd build
for GOOS in darwin linux windows; do
    for GOARCH in 386 amd64; do
        echo "`GOOS=${GOOS} GOARCH=${GOARCH} go build -gcflags="-l=4" -compiler="gc" -v -o nucamino-${GOOS}-${GOARCH} github.com/hivdb/nucamino/cmd/nucamino 2>&1` => ./build/nucamino-${GOOS}-${GOARCH}"
    done
done
mv nucamino-windows-amd64 nucamino-windows-amd64.exe
mv nucamino-windows-386 nucamino-windows-386.exe
