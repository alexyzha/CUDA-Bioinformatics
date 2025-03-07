all: build run

build:
	docker build --platform linux/amd64 -t cuda-qcb .

run:
	docker run --rm --gpus all --platform linux/amd64 -it -v $(shell pwd):/src cuda-qcb

prune:
	docker container prune -f

delete:
	docker rmi -f cuda-qcb
