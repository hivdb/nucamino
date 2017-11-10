#! /bin/sh

set -e
/go/src/github.com/hivdb/nucamino/hack/build.sh _pprof
cd /go/src/github.com/hivdb/nucamino
./_pprof/nucamino hiv1b --pprof --gene POL --input _pprof_input.fas -o /dev/null 2> /tmp/output
PPROF_FILE=`grep "\(cpu\|memory\) profiling disabled" /tmp/output | awk '{print $(NF)}'`
mkdir -p _local
if [ "$1" = "pdf" ]; then
    go tool pprof -pdf ${PPROF_FILE} > _local/pprof.pdf
else
    go tool pprof ${PPROF_FILE}
fi
echo "Profile result generated at ./_local/pprof.pdf"
# clean up
rm -rf _pprof
