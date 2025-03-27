$(shell mkdir -p $(./obj))

target: clean cpu gpu

cpu:
	g++ -std=c++17 helloworld.cpp -o ./obj/helloworld

gpu:
	nvcc helloworld.cu -o ./obj/cu_helloworld

clean:
	rm -rf $(./obj)/*.o $(helloworld cu_helloworld)
