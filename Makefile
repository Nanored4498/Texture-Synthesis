EXEC=main

FLAGS=-Wall -O3 -fopenmp

BUILD_FOLDER=build
SRC_C=$(wildcard *.c)
OBJ_C=$(patsubst %.c, $(BUILD_FOLDER)/%.o, $(SRC_C))
SRC_CPP=$(wildcard *.cpp)
OBJ_CPP=$(patsubst %.cpp, $(BUILD_FOLDER)/%.opp, $(SRC_CPP))

all: $(BUILD_FOLDER) $(EXEC)

$(BUILD_FOLDER):
	mkdir $(BUILD_FOLDER)

$(EXEC): $(OBJ_C) $(OBJ_CPP)
	g++ $(FLAGS) $^ -o main

$(BUILD_FOLDER)/%.o: %.c
	gcc $(FLAGS) -c $< -o $@

$(BUILD_FOLDER)/%.opp: %.cpp
	g++ $(FLAGS) -c $< -o $@