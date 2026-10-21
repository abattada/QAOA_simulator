# Compiler and flags
CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra

# Layout
SRC_DIR   = src
TEST_DIR  = tests
BUILD_DIR = build
SOURCES   = $(wildcard $(SRC_DIR)/*.cpp)
HEADERS   = $(wildcard $(SRC_DIR)/*.hpp)
OBJECTS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
LIB_OBJS  = $(filter-out $(BUILD_DIR)/main.o,$(OBJECTS))

# Targets
TARGET    = qaoa_main
TEST_BIN  = $(BUILD_DIR)/test_qaoa

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TEST_BIN): $(TEST_DIR)/test_qaoa.cpp $(LIB_OBJS) $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -I$(SRC_DIR) $< $(LIB_OBJS) -o $@

# Run from the repo root: the tests read graph/qml_bi_4
test: $(TEST_BIN)
	./$(TEST_BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all test clean
