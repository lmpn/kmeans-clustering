SHELL = /bin/sh
CXX= g++
LD  =g++
CXXFLAGS= -Xpreprocessor -fopenmp -O3 -std=c++11 -ggdb3 -Wno-unused-parameter

ifneq ($(MACOS), yes)
	LIBS = -I/share/apps/papi/5.5.0/include -L/share/apps/papi/5.5.0/lib
	CXXFLAGS+= -lpapi
endif

 

BIN_NAME = kmeans
SRC_DIR = src
BIN_DIR = bin
HEADER_DIR = header
BUILD_DIR = build
ASM_DIR = asm
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(patsubst src/%.cpp,build/%.o,$(SRC))
ASM = $(patsubst src/%.cpp,asm/%.s,$(SRC))
INCLUDES = -I$(HEADER_DIR)
vpath %.cpp $(SRC_DIR)
#######
#RULES
#######
.DEFAULT_GOAL = all


$(BIN_DIR)/$(BIN_NAME):  $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJ)
$(BUILD_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@ 
	
checkdirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

all: checkdirs $(BIN_DIR)/$(BIN_NAME)

clean:
	rm -rf $(BUILD_DIR)/* $(BIN_DIR)/* 
