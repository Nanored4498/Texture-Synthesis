EXEC=main

FLAGS=-Wall -O3 -fopenmp

BUILD_FOLDER=build
SRC_C=$(wildcard *.c)
OBJ_C=$(patsubst %.c, $(BUILD_FOLDER)/%.o, $(SRC_C))
SRC_CPP=$(wildcard *.cpp)
OBJ_CPP=$(patsubst %.cpp, $(BUILD_FOLDER)/%.opp, $(SRC_CPP))

DEBUG_FOLDER=debug
DEB=$(DEBUG_FOLDER)/main
DEB_C=$(patsubst %.c, $(DEBUG_FOLDER)/%.o, $(SRC_C))
DEB_CPP=$(patsubst %.cpp, $(DEBUG_FOLDER)/%.opp, $(SRC_CPP))

all: $(BUILD_FOLDER) $(EXEC)

$(BUILD_FOLDER):
	mkdir $(BUILD_FOLDER)

$(EXEC): $(OBJ_C) $(OBJ_CPP)
	g++ $(FLAGS) $^ -o $@

$(BUILD_FOLDER)/%.o: %.c
	gcc $(FLAGS) -c $< -o $@

$(BUILD_FOLDER)/%.opp: %.cpp
	g++ $(FLAGS) -c $< -o $@


#DEBUG PART

.PHONY: deb
deb: $(DEBUG_FOLDER) $(DEB)
	valgrind ./$(DEB) ${ARGS}

$(DEBUG_FOLDER):
	mkdir $(DEBUG_FOLDER)

$(DEB): $(DEB_C) $(DEB_CPP)
	g++ $(FLAGS) -g $^ -o $@

$(DEBUG_FOLDER)/%.o: %.c
	gcc $(FLAGS) -g -c $< -o $@

$(DEBUG_FOLDER)/%.opp: %.cpp
	g++ $(FLAGS) -g -c $< -o $@


#CLEAN PART

clean:
	rm -r $(BUILD_FOLDER) $(DEBUG_FOLDER)