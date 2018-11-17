SHELL = /bin/sh
CXX= g++
LD  =g++
CXXFLAGS= -Xpreprocessor -fopenmp -O3 -std=c++11 -g -ggdb3 -Wno-unused-parameter
ifeq($(VERSION), seq)
	BIN_NAME = kmeans_si
endif

ifeq($(VERSION), seqp)
	LIBS = -I/share/apps/papi/5.5.0/include -L/share/apps/papi/5.5.0/lib
	CXXFLAGS+= -lpapi -DPAPI
	BIN_NAME = kmeans_spi
endif

ifeq($(VERSION), seqfp)
	LIBS = -I/share/apps/papi/5.5.0/include -L/share/apps/papi/5.5.0/lib
	CXXFLAGS+= -lpapi -DPAPI
	BIN_NAME = kmeans_spf
endif

ifeq($(VERSION), seqfs)
	BIN_NAME = kmeans_sfs
endif

ifeq($(VERSION), seqf)
	BIN_NAME = kmeans_sf
endif
ifeq ($(MACOS), no)
	LIBS = -I/share/apps/papi/5.5.0/include -L/share/apps/papi/5.5.0/lib
	CXXFLAGS+= -lpapi -DPAPI
endif

BIN_NAME = kmeans_spi
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj
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


$(BIN_DIR)/$(BIN_NAME):  $(OBJ) $(ASM)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -o $@ $(OBJ)
$(BUILD_DIR)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@ 
$(ASM_DIR)/%.s: %.cpp
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $(LIBS) $< -o $@ 
	
checkdirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(ASM_DIR)

all: checkdirs $(BIN_DIR)/$(BIN_NAME)

delete:
	rm -rf $(BUILD_DIR) $(BIN_DIR) $(ASM_DIR) $(OBJ_DIR)
clean:
	rm -rf $(BUILD_DIR)/* $(BIN_DIR)/* $(ASM_DIR)/* $(OBJ_DIR)/* out outp
