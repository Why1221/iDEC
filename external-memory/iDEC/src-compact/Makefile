CXX = g++
CXXFLAGS = -std=c++11 -O3 -w -I../../Common
LDFLAGS =
RM = gio trash -f 

TARGETS = idec idec-steal rtree-example


INT_SRC = $(wildcard *.cpp)
INT_OBJ = $(INT_SRC:.cpp=.o)

TARGET_SRC = $(wildcard *.cc)
TAR_OBJ = $(TARGET_SRC:.cc=.o)

DIFF_SRC = $(wildcard *.cxx)
DIFF_OBJ = $(DIFF_SRC:.cxx=.o)


all: $(TARGETS)


rtree-example: rtree-example.o $(INT_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

idec: main.o $(INT_OBJ) $(DIFF_OBJ)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

idec-steal: main.cc $(INT_SRC) $(DIFF_SRC)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) -DUSE_SEARCH_HISTORY=1


.PHONY: clean
clean:
	$(RM) $(INT_OBJ) $(TAR_OBJ) $(TARGETS) $(DIFF_OBJ)

