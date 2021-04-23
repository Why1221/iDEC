CXX = g++
# please disable all warings if you want to use un-copmpact version 
# CXXFLAGS = -std=c++11 -O3 -w
CXXFLAGS = -std=c++11 -O3 -w -I../Common -DMEASURE_IO_TIME
LDFLAGS = -lgflags
RM = gio trash -f

TARGETS = rtree-example idec-compact idec idec-steal


INT_SRC_COMP = $(wildcard src-compact/*.cpp)
INT_OBJ_COMP = $(INT_SRC_COMP:.cpp=.o)

DIFF_SRC_COMP = $(wildcard src-compact/*.cxx)
DIFF_OBJ_COMP = $(DIFF_SRC_COMP:.cxx=.o)

TARGET_SRC_COMP = $(wildcard src-compact/*.cc)
TAR_OBJ_COMP = $(TARGET_SRC_COMP:.cc=.o)

INT_SRC = $(wildcard src/*.cpp)
INT_OBJ = $(filter-out src/main-compact.o src/srs-compact.o, $(INT_SRC:.cpp=.o))
TARGET_SRC = src/main.cc
TAR_OBJ = $(TARGET_SRC:.cc=.o)



all: $(TARGETS)

# TODO: figure out why cxx file cannot work 
src-compact/idec.o: $(DIFF_SRC_COMP)
	-$(CXX) -c -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

rtree-example: src-compact/rtree-example.o $(INT_OBJ_COMP)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

idec-compact: src-compact/main.o $(INT_OBJ_COMP) $(DIFF_OBJ_COMP)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

idec-steal: src-compact/main.cc $(INT_SRC_COMP) $(DIFF_SRC_COMP)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -DUSE_SEARCH_HISTORY

idec: src/main.o $(INT_OBJ)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)


.PHONY: clean
clean:
	-$(RM) $(INT_OBJ_COMP) $(TAR_OBJ_COMP) $(INT_OBJ) $(TAR_OBJ) $(TARGETS) $(DIFF_OBJ_COMP)