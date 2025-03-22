## Ubuntu base image
FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

## Install basic
RUN apt-get update && apt-get install -y \
    sra-toolkit \
    build-essential \
    valgrind \
    gdb \
    cmake \
    wget \
    curl \
    git \
    unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

## Install GTEST
RUN git clone -q https://github.com/google/googletest.git /googletest \
    && mkdir -p /googletest/build \
    && cd /googletest/build \
    && cmake .. && make -j$(nproc) && make install \
    && ldconfig \
    && rm -rf /googletest

## SETWD
WORKDIR /src/

## Bash
CMD ["/bin/bash"]
