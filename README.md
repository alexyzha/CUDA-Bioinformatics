# About This Project

All documentation can be found [here](https://github.com/alexyzha/CUDA-Bioinformatics/wiki). 

## Why?

I decided to write a bioinformatics toolchain because even though I am a Quantitative/Computational Biology major, we don't get to work with that much code in class. The things I've done with code so far have mostly been statistics/data related, and I have a lot more fun coding from the ground up. Additionally, I wanted to get more experience coding in `CUDA`, and algorithms used in bioinformatics are 1. massively parallel and 2. somewhat familiar to me. Finally, I also wanted to gain more experience setting up things like unit/integration tests and CI/CD pipelines. 

## What's in it? 

The code in this repository is written to parse, format/package, and analyze `.fasta`, `.fastq`, and `.sam` files. The regular `C++` code supports all 3 of these files. However, the `CUDA` code only supports operations on `.fastq` file data. However, all of the `CUDA` code is neatly wrapped in `__host__` code wrappers: you don't actually have to process/format data to use the `CUDA` kernels I wrote.

## Developing Environment

- Docker base image: `cuda:12.1.1-devel-ubuntu22.04`. On my Mac, I use `Ubuntu:LATEST` since I don't have an NVIDIA GPU on there.
- Docker packages: `GTEST`, `sra-toolkit`, `valgrind`, `gdb`, `cmake`, `build-essential`, and some other not super important ones (see `./Dockerfile`).
- On my desktop, I use `nvidia-ctk` to run `CUDA` code inside a Docker container with a `RTX 4060`.