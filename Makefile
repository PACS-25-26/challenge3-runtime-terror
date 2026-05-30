CXX = mpicxx # Use mpicxx for MPI support
CXXFLAGS = -Wall -Wextra -O3 -std=c++17 -fopenmp -DOMPI_SKIP_MPICXX=1
INCLUDES = -I./include

SRC_DIR = src
OBJ_DIR = obj

# Find all .cpp files in the src directory
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
# Generate corresponding .o files in the obj directory
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

TARGET = main.out

all: directories $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

directories:
	@mkdir -p $(OBJ_DIR)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean directories