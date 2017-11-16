#!/bin/sh
# modified based on: https://github.com/codecov/example-go

set -e
echo "" > coverage.txt
cd `dirname $0`/..

for d in $(go list ./... | \grep -v vendor); do
    if [ -f /etc/alpine-release ]; then
        go test -coverprofile=profile.out -covermode=atomic $d
    else
        go test -race -coverprofile=profile.out -covermode=atomic $d
    fi
    if [ -f profile.out ]; then
        cat profile.out >> coverage.txt
        rm profile.out
    fi
done
