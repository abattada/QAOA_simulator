# Compiler and flags
CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra

# Layout
SRC_DIR   = src
BUILD_DIR = build
SOURCES   = $(wildcard $(SRC_DIR)/*.cpp)
HEADERS   = $(wildcard $(SRC_DIR)/*.hpp)
OBJECTS   = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

# Target executable
TARGET = qaoa_main

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(HEADERS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR) $(TARGET)

.PHONY: all clean
