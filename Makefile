APPROOT = /go/src/github.com/hivdb/nucamino
GOOS ?= $(shell uname | tr A-Z a-z)
GOARCH ?= $(subst x86_64,amd64,$(patsubst i%86,386,$(shell uname -m)))


dockerimage:
	@docker build . -qt nucamino-dev

pprof: dockerimage
	@docker rm -f nucamino-pprof 2>/dev/null || true
	@docker run --rm -it --name nucamino-pprof --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/hack/pprof.sh

pprof-pdf: dockerimage
	@docker rm -f nucamino-pprof 2>/dev/null || true
	@docker run --rm -it --name nucamino-pprof --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/hack/pprof.sh pdf

build: dockerimage
	@docker rm -f nucamino-build 2>/dev/null || true
	@docker run --rm -it --name nucamino-build --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} -e GOOS=${GOOS} -e GOARCH=${GOARCH} nucamino-dev ${APPROOT}/hack/build.sh

cross-build: dockerimage
	@docker rm -f nucamino-build 2>/dev/null || true
	@docker run --rm -it --name nucamino-build --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/hack/cross-compile.sh

.PHONY: dockerimage pprof
