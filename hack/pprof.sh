#! /bin/sh

set -e

if [ -z $GOPATH ]; then
    GOPATH=${HOME}/go
fi

${GOPATH}/src/github.com/hivdb/nucamino/hack/build.sh _pprof
cd ${GOPATH}/src/github.com/hivdb/nucamino
./_pprof/nucamino align hiv1b POL --pprof -i _pprof_input.fas -o /dev/null 2> /tmp/output
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
