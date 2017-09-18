#! /bin/sh

set -e

cd /go/src/github.com/hivdb/nucamino
go get ./...
cd pprof
go build main.go
./main 2> /tmp/output
PPROF_FILE=`grep "\(cpu\|memory\) profiling disabled" /tmp/output | awk '{print $(NF)}'`
mkdir -p ../local
if [ "$1" = "pdf" ]; then
    go tool pprof -runtime -pdf ${PPROF_FILE} > ../local/pprof.pdf
else
    go tool pprof -runtime ${PPROF_FILE}
fi
echo "Profile result generated at ./local/pprof.pdf"
# clean up
rm main
