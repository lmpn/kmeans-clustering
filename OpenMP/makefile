CC= g++-6
CFLAGS= -Xpreprocessor -fopenmp -lomp -O3 -std=c++11  -lm -g -Wno-unused-parameter
SRC_DIR = src/
BIN_DIR = bin/
HEADER_DIR = header/
OBJ_DIR = obj/
SRC = $(wildcard $(SRC_DIR)*.cpp)
OBJ = $(patsubst src/%.cpp,obj/%.o,$(SRC))

ASM = $(patsubst src/%.cpp,asm/%.s,$(SRC))
.DEFAULT_GOAL = all



all: checkdirs $(OBJ) $(ASM) obj/main.o asm/main.s
	$(CC) -o $(BIN_DIR)run -I$(HEADER_DIR) $(CFLAGS) obj/main.o

obj/main.o: main.cpp
	$(CC) -c -I$(HEADER_DIR) $(CFLAGS) $< -o $@

asm/main.s: main.cpp
	$(CC) -S -I$(HEADER_DIR) $(CFLAGS) $< -o $@

#Create Assembly
$(ASM): $(SRC)
	$(CC) -S -I$(HEADER_DIR) $(CFLAGS) $< -o $@

#Create Objects
$(OBJ): $(SRC)
	$(CC) -c -I$(HEADER_DIR) $(CFLAGS) $< -o $@


checkdirs:
	@mkdir -p build
	@mkdir -p bin

install:
	mkdir -p src bin obj asm header

.PHONY: clean

clean:
	rm -rf bin/* obj/* asm/*
