FROM ubuntu:20.04

MAINTAINER Gregor Rot <gregor.rot@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Zurich

RUN apt update
RUN echo '8'| apt-get install -y tzdata
RUN apt-get install -y python3-pip xxd wget samtools python3-dev zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libcurl4-openssl-dev libssl-dev
RUN pip install pybio

WORKDIR /
RUN wget https://github.com/alexdobin/STAR/archive/2.7.11a.tar.gz
RUN tar -xzf 2.7.11a.tar.gz
WORKDIR /STAR-2.7.11a/source
RUN make

ENV PATH=$PATH:/STAR-2.7.11a/source
