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
	@open _local/pprof.pdf

build: dockerimage
	@docker rm -f nucamino-build 2>/dev/null || true
	@docker run --rm -it --name nucamino-build --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} -e GOOS=${GOOS} -e GOARCH=${GOARCH} nucamino-dev ${APPROOT}/hack/build.sh

cross-build: dockerimage
	@docker rm -f nucamino-build 2>/dev/null || true
	@docker run --rm -it --name nucamino-build --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/hack/cross-compile.sh

shell: dockerimage
	@docker rm -f nucamino-build 2>/dev/null || true
	@docker run --rm -it --name nucamino-build --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev /bin/sh

test: dockerimage
	@docker rm -f nucamino-test 2>/dev/null || true
	@docker run --rm -it --name nucamino-test --volume $(shell pwd)/.goenv:/go --volume $(shell pwd):${APPROOT} nucamino-dev ${APPROOT}/hack/test.sh

aws_lambda_zip: cross-build
	@mkdir -p /tmp/_nucamino_aws_lambda_zip
	@cp build/nucamino-linux-amd64 /tmp/_nucamino_aws_lambda_zip/nucamino
	@chmod +x /tmp/_nucamino_aws_lambda_zip/nucamino
	@cp scripts/aws_lambda_handler.py /tmp/_nucamino_aws_lambda_zip/nucamino_handler.py
	@cd /tmp/_nucamino_aws_lambda_zip && zip nucamino-aws-lambda.zip nucamino nucamino_handler.py
	@mv /tmp/_nucamino_aws_lambda_zip/nucamino-aws-lambda.zip build/
	@ls build/nucamino-aws-lambda.zip
	@rm -rf /tmp/_nucamino_aws_lambda_zip

.PHONY: dockerimage pprof
