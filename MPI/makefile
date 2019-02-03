SHELL		:= /bin/sh
CXX 		:= mpicxx
LD			:= mpicxx
CXXFLAGS 	:= -O3 -std=c++11 -g -ggdb3 -Wno-unused-parameter -DREDUCEBCAST
BIN_NAME 	:= rbcore

LIB			:= -Llib
BIN 		:= bin
BUILD 		:= build
INCLUDES 	:= -Iinclude
SRC 		:= $(wildcard src/*.cpp)
OBJ 		:= $(patsubst src/%.cpp,build/%.o,src/)
ASM 		:= $(patsubst src/%.cpp,asm/%.s,src/)
vpath %.cpp $(SRC)




#######
#RULES
#######
.DEFAULT_GOAL = all


$(BIN)/$(BIN_NAME):  $(OBJ) $(ASM)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJ)
$(BUILD)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@ 
$(ASM)/%.s: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@ 
checkdirs:
	@mkdir -p $(BUILD)
	@mkdir -p $(BIN)
	@mkdir -p $(OBJ)
	@mkdir -p $(ASM
all: checkdirs $(BIN_DIR)/$(BIN_NAME)

delete:
	rm -rf $(BUILD_DIR) $(BIN_DIR) $(ASM_DIR) $(OBJ_DIR)
clean:
	rm -rf $(BUILD_DIR)/* $(ASM_DIR)/* $(OBJ_DIR)/* 
