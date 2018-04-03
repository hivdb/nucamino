FROM golang:1.10-stretch
RUN apt-get update -q && apt-get install -qy git graphviz fonts-noto
VOLUME ["/go"]
