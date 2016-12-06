APPROOT = /go/src/github.com/hivdb/nucamino

dockerimage:
	@docker build . -qt nucamino-dev

pprof: dockerimage
	@docker rm -f nucamino-pprof 2>/dev/null || true
	@docker run --rm -it --name nucamino-pprof --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/pprof/pprof.sh

pprof-pdf: dockerimage
	@docker rm -f nucamino-pprof 2>/dev/null || true
	@docker run --rm -it --name nucamino-pprof --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/pprof/pprof.sh pdf

.PHONY: dockerimage pprof
