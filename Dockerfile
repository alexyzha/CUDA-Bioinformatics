## Ubuntu base image
FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive

## Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    cmake \
    wget \
    curl \
    git \
    unzip \
    && rm -rf /var/lib/apt/lists/*

## Install GTEST
RUN git clone -q https://github.com/google/googletest.git /googletest \
    && mkdir -p /googletest/build \
    && cd /googletest/build \
    && cmake .. && make && make install \
    && cd / && rm -rf /googletest

## Setwd
WORKDIR /src/

## Bash
CMD ["/bin/bash"]
