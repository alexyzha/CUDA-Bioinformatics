########################################
##        ALL CPU MODULE TESTS        ##
##   COMPILES [./obj/all_cpu_tests]   ##
##         CAN AUTO RUN IN GA         ##
########################################

CXX 	 	= g++
CXXFLAGS 	= -std=c++17 -I../../src -I../../src/cpu -Wno-trigraphs
WARNINGS 	= -Wall -Wextra
LDFLAGS  	= -lgtest -lgtest_main -pthread
SRC_DIR  	= ../../src
CPU_DIR  	= ../../src/cpu
OBJ_DIR  	= ./obj
TEST_DIR 	= .

$(shell mkdir -p $(OBJ_DIR))

SRC_FILES = \
	$(CPU_DIR)/fa_read.cpp \
	$(CPU_DIR)/fq_read.cpp \
	$(CPU_DIR)/fx_parser.cpp \
	$(CPU_DIR)/fx_util.cpp \
	$(CPU_DIR)/sam_container.cpp \
	$(CPU_DIR)/util_structs.cpp \
	$(CPU_DIR)/util.cpp \
	$(CPU_DIR)/xam_parser.cpp \
	$(CPU_DIR)/xam_util.cpp

TEST_FILES = \
	$(TEST_DIR)/fa_read_tests.cpp \
	$(TEST_DIR)/fq_read_tests.cpp \
	$(TEST_DIR)/fx_parser_tests.cpp \
	$(TEST_DIR)/fx_util_tests.cpp \
	$(TEST_DIR)/sam_container_tests.cpp \
	$(TEST_DIR)/sam_parser_tests.cpp \
	$(TEST_DIR)/sam_read_tests.cpp \
	$(TEST_DIR)/util_structs_tests.cpp \
	$(TEST_DIR)/util_tests.cpp \
	$(TEST_DIR)/xam_util_tests.cpp

OBJ_FILES = $(patsubst %.cpp, $(OBJ_DIR)/%.o, $(notdir $(SRC_FILES) $(TEST_FILES)))
OUTPUT = $(OBJ_DIR)/all_cpu_tests

target: clean test

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
