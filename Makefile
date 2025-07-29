BIN_DIR := bin
BUILD_DIR := build
SRC_DIRS := src protocol/src utils/src

CXX := /opt/homebrew/opt/llvm/bin/clang++
CXXFLAGS := -Wall -march=native -O3 -std=c++17 -MMD -MP
INCLUDES := -I/opt/homebrew/opt/llvm/include \
            -I/opt/homebrew/include/eigen3 \
            -I/opt/homebrew/include \
            -I.

LDFLAGS := -L/opt/homebrew/opt/llvm/lib \
           -L/opt/homebrew/lib
LDLIBS := -lhdf5 -lhdf5_cpp

SOURCES := $(wildcard src/*.cpp) $(wildcard protocol/src/*.cpp) $(wildcard utils/src/*.cpp)
OBJECTS := $(patsubst %.cpp, $(BUILD_DIR)/%.o, $(notdir $(SOURCES)))

.PHONY: all
all: builddirs $(OBJECTS)

.PHONY: builddirs
builddirs:
	mkdir -p $(BIN_DIR)
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: src/%.cpp | builddirs
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
$(BUILD_DIR)/%.o: protocol/src/%.cpp | builddirs
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
$(BUILD_DIR)/%.o: utils/src/%.cpp | builddirs
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

.PHONY: build
build: all
	@echo "Run \`make run EXE=<file.cpp>\` to compile a specific main file."

.PHONY: run
run: build
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $(BIN_DIR)/$(basename $(notdir $(EXE))) main/$(EXE) $(OBJECTS) $(LDLIBS)

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

-include $(OBJECTS:.o=.d)
