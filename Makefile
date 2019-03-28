EXEC=main
GUI=gui

FLAGS=-Wall -O3 -fopenmp
LINK=-lboost_filesystem -lboost_system
LINK_GTK=`pkg-config gtkmm-3.0 --cflags --libs`

BUILD_FOLDER=build
SRC_C=$(wildcard *.c)
OBJ_C=$(patsubst %.c, $(BUILD_FOLDER)/%.o, $(SRC_C))
SRC_CPP=$(wildcard *.cpp)
OBJ_CPP=$(patsubst %.cpp, $(BUILD_FOLDER)/%.opp, $(SRC_CPP))
OBJ_CPP_GUI=$(filter-out $(BUILD_FOLDER)/main.opp, $(OBJ_CPP)) $(BUILD_FOLDER)/gui.opp

DEBUG_FOLDER=debug
DEB=$(DEBUG_FOLDER)/main
DEB_C=$(patsubst %.c, $(DEBUG_FOLDER)/%.o, $(SRC_C))
DEB_CPP=$(patsubst %.cpp, $(DEBUG_FOLDER)/%.opp, $(SRC_CPP))

IMS=magnific.png map_*.png out_*.png

all: $(BUILD_FOLDER) $(EXEC)

g: $(BUILD_FOLDER) $(GUI)

$(BUILD_FOLDER):
	mkdir $(BUILD_FOLDER)

$(EXEC): $(OBJ_C) $(OBJ_CPP)
	g++ $(FLAGS) $^ -o $@ $(LINK)

$(GUI): $(OBJ_C) $(OBJ_CPP_GUI)
	g++ $(FLAGS) $^ -o $@ $(LINK) $(LINK_GTK)

$(BUILD_FOLDER)/%.o: %.c
	gcc $(FLAGS) -c $< -o $@

$(BUILD_FOLDER)/gui.opp: GUI/gui.cpp
	g++ $(FLAGS) -c $< -o $@ $(LINK_GTK)

$(BUILD_FOLDER)/%.opp: %.cpp
	g++ $(FLAGS) -c $< -o $@


#DEBUG PART

.PHONY: deb
deb: $(DEBUG_FOLDER) $(DEB)
	valgrind ./$(DEB) ${ARGS}

$(DEBUG_FOLDER):
	mkdir $(DEBUG_FOLDER)

$(DEB): $(DEB_C) $(DEB_CPP)
	g++ $(FLAGS) -g $^ -o $@ $(LINK)

$(DEBUG_FOLDER)/%.o: %.c
	gcc $(FLAGS) -g -c $< -o $@

$(DEBUG_FOLDER)/%.opp: %.cpp
	g++ $(FLAGS) -g -c $< -o $@


#CLEAN PART

clean:
	rm -rf $(BUILD_FOLDER) $(DEBUG_FOLDER)
clean_im:
	rm -f $(IMS)
mr_proper: clean clean_im
	rm -f $(GUI) $(EXEC)