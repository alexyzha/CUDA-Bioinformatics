########################################
##        ALL GPU MODULE TESTS        ##
##   COMPILES [./obj/all_gpu_tests]   ##
##         CAN AUTO RUN IN GA         ##
########################################

NVCC        = nvcc
CXX         = g++
CXXFLAGS    = -std=c++17 -I../../src -I../../src/cpu -I../../src/gpu -I../../src/gpu/headers -I./headers
NVCCFLAGS   = -std=c++17 -dc -rdc=true -I../../src -I../../src/cpu -I../../src/gpu -I../../src/gpu/headers -I./headers
OBJ_DIR     = ./obj
SRC_DIR     = ../../src
CPU_DIR     = $(SRC_DIR)/cpu
GPU_DIR     = $(SRC_DIR)/gpu
TEST_DIR    = .
OUTPUT      = $(OBJ_DIR)/all_gpu_tests

$(shell mkdir -p $(OBJ_DIR))

# CPU core
SRC_FILES = \
	$(CPU_DIR)/fa_read.cpp \
	$(CPU_DIR)/fq_read.cpp \
	$(CPU_DIR)/fx_parser.cpp \
	$(CPU_DIR)/fx_util.cpp \
	$(CPU_DIR)/sam_container.cpp \
	$(CPU_DIR)/util_structs.cpp \
	$(CPU_DIR)/util.cpp \
	$(CPU_DIR)/xam_parser.cpp \
	$(CPU_DIR)/xam_util.cpp \
	$(TEST_DIR)/test_util.cpp

# GPU core
GPU_SRC_FILES = \
	$(GPU_DIR)/fq_filters.cu \
	$(GPU_DIR)/kmer_util.cu \
	$(GPU_DIR)/local_alignment.cu \
	$(GPU_DIR)/util_structs.cu \
	$(GPU_DIR)/util.cu \
	$(GPU_DIR)/wrappers.cu

# GPU test files
GPU_TEST_FILES = \
	$(TEST_DIR)/all_gpu_tests.cu \
	$(TEST_DIR)/global_wrappers.cu \
	$(TEST_DIR)/util_structs_tests.cu \
	$(TEST_DIR)/util_tests.cu \
	$(TEST_DIR)/wrapper_tests.cu

CPP_OBJS = $(patsubst $(CPU_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(filter $(CPU_DIR)/%, $(SRC_FILES)))
CPP_OBJS += $(patsubst $(TEST_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(filter $(TEST_DIR)/%, $(SRC_FILES)))
CU_OBJS  = $(patsubst $(GPU_DIR)/%.cu, $(OBJ_DIR)/%.o, $(GPU_SRC_FILES))
CU_OBJS += $(patsubst $(TEST_DIR)/%.cu, $(OBJ_DIR)/%.o, $(GPU_TEST_FILES))
OBJ_FILES = $(CPP_OBJS) $(CU_OBJS)

all: $(OUTPUT)

# Compile .cpp files
$(OBJ_DIR)/%.o: $(CPU_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cu files
$(OBJ_DIR)/%.o: $(GPU_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Link all
$(OUTPUT): $(OBJ_FILES)
	$(NVCC) -rdc=true $(CXXFLAGS) $^ -o $@

test: all
	./$(OUTPUT)

clean:
	rm -rf $(OBJ_DIR)/*.o $(OUTPUT)