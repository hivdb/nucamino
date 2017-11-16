FROM golang:1.9-alpine3.6
RUN apk add --no-cache git graphviz ttf-droid gcc
VOLUME ["/go"]
