EXEC=main
GUI=gui

FLAGS=-Wall -O3 -fopenmp
LINK=-lboost_filesystem -lboost_system
LINK_GTK=`pkg-config gtkmm-3.0 --cflags --libs`
INCLUDE=-I src/headers

BUILD_FOLDER=build
SRC_C=$(wildcard src/*.c)
OBJ_C=$(patsubst src/%.c, $(BUILD_FOLDER)/%.o, $(SRC_C))
SRC_CPP=$(wildcard src/*.cpp)
OBJ_CPP=$(patsubst src/%.cpp, $(BUILD_FOLDER)/%.opp, $(SRC_CPP))
OBJ_CPP_GUI=$(filter-out $(BUILD_FOLDER)/main.opp, $(OBJ_CPP)) $(BUILD_FOLDER)/gui.opp

DEBUG_FOLDER=debug
DEB=$(DEBUG_FOLDER)/main
DEB_C=$(patsubst src/%.c, $(DEBUG_FOLDER)/%.o, $(SRC_C))
DEB_CPP=$(patsubst src/%.cpp, $(DEBUG_FOLDER)/%.opp, $(SRC_CPP))

IMS=magnific.png map_*.png out*.png

all: $(BUILD_FOLDER) $(EXEC)

g: $(BUILD_FOLDER) $(GUI)

$(BUILD_FOLDER):
	mkdir $(BUILD_FOLDER)

$(EXEC): $(OBJ_C) $(OBJ_CPP)
	g++ $(FLAGS) $^ -o $@ $(LINK)

$(GUI): $(OBJ_C) $(OBJ_CPP_GUI)
	g++ $(FLAGS) $^ -o $@ $(LINK) $(LINK_GTK)

$(BUILD_FOLDER)/%.o: src/%.c
	gcc $(FLAGS) -c $< -o $@ $(INCLUDE)

$(BUILD_FOLDER)/gui.opp: src/GUI/gui.cpp
	g++ $(FLAGS) -c $< -o $@ $(LINK_GTK) $(INCLUDE)

$(BUILD_FOLDER)/%.opp: src/%.cpp
	g++ $(FLAGS) -c $< -o $@ $(INCLUDE)


#DEBUG PART

.PHONY: deb
deb: $(DEBUG_FOLDER) $(DEB)
	valgrind ./$(DEB) ${ARGS}

$(DEBUG_FOLDER):
	mkdir $(DEBUG_FOLDER)

$(DEB): $(DEB_C) $(DEB_CPP)
	g++ $(FLAGS) -g $^ -o $@ $(LINK)

$(DEBUG_FOLDER)/%.o: src/%.c
	gcc $(FLAGS) -g -c $< -o $@ $(INCLUDE)

$(DEBUG_FOLDER)/%.opp: src/%.cpp
	g++ $(FLAGS) -g -c $< -o $@ $(INCLUDE)


#CLEAN PART

clean:
	rm -rf $(BUILD_FOLDER) $(DEBUG_FOLDER)
clean_im:
	rm -f $(IMS)
mr_proper: clean clean_im
	rm -f $(GUI) $(EXEC)