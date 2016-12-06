FROM golang:1.7.4-alpine
RUN apk add --no-cache git graphviz ttf-droid
VOLUME ["/go"]
