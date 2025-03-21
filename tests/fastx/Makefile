########################################
##         FASTX MODULE TESTS         ##
##     COMPILES [./obj/all_tests]     ##
##         CAN AUTO RUN IN GA         ##
########################################

CXX = g++
CXXFLAGS = -std=c++17 -I../../src -I../../src/cpu
WARNINGS = -Wall -Wextra
LDFLAGS = -lgtest -lgtest_main -pthread
SRC_DIR = ../../src
CPU_DIR = ../../src/cpu
OBJ_DIR = ./obj
TEST_DIR = .

$(shell mkdir -p $(OBJ_DIR))

SRC_FILES = $(SRC_DIR)/fa_read.cpp $(SRC_DIR)/fq_read.cpp $(SRC_DIR)/util.cpp $(CPU_DIR)/fx_parser.cpp $(SRC_DIR)/util.cpp
TEST_FILES = $(TEST_DIR)/fa_read_tests.cpp $(TEST_DIR)/fq_read_tests.cpp $(TEST_DIR)/fx_parser_tests.cpp
OBJ_FILES = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES) $(TEST_FILES)))
OUTPUT = $(OBJ_DIR)/all_tests

all: $(OUTPUT)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(CPU_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OUTPUT): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

test: all
	./$(OUTPUT)

clean:
	rm -rf $(OBJ_DIR)/*.o $(OUTPUT)
