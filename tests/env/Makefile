OBJ_DIR = ./obj

$(shell mkdir -p $(OBJ_DIR))

target: clean cpu gpu

cpu:
	g++ -std=c++17 helloworld.cpp -o $(OBJ_DIR)/helloworld

gpu:
	nvcc helloworld.cu -o $(OBJ_DIR)/cu_helloworld

clean:
	rm -rf $(OBJ_DIR)/*.o $(OBJ_DIR)/helloworld $(OBJ_DIR)/cu_helloworld
