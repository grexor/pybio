FROM ubuntu:20.04

MAINTAINER Gregor Rot <gregor.rot@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

RUN apt update
RUN echo '8'| apt-get install -y tzdata
RUN apt-get install -y python3-pip wget samtools rna-star python3-dev zlib1g-dev libbz2-dev libl
zma-dev libncurses5-dev libcurl4-openssl-dev libssl-dev
RUN pip install pybio
