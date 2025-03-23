all: run

run:
	docker run --runtime=nvidia --gpus all --platform linux/amd64 -it --rm -v $(shell pwd)/../:/src cuda-qcb
