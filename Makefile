all: build run

build:
	docker build --platform linux/amd64 -t cuda-qcb .

run:
	docker run --platform linux/amd64 -it --rm -v $(PWD):/src cuda-qcb

prune:
	docker container prune -f

delete:
	docker rmi -f cuda-qcb

linux:
	docker run --runtime=nvidia --gpus all --platform linux/amd64 -it --rm -v $(shell pwd)/:/src cuda-qcb
