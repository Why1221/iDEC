CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -I ./src -I ../Common 
# CXXFLAGS = -O3 -std=c++11 -Wall -I ./src -I ../Common -DDEBUG=1 -DSAVE_PROJ_DATA=1
LDFALAGS = 
LDFALAGS_TEST = -lgtest -lpthread
RM = gio trash -f 

TARGETS = idec idec-steal
SRCS = $(wildcard src/*.cpp)
OBJS = $(SRCS:.cpp=.o)
TAR_SRCS = $(wildcard src/*.cc)
TAR_OBJS = $(TAR_SRCS:.cc=.o)

# for tests
TEST_SRCS = $(wildcard test/*.cc)
TEST_OBJS = $(TEST_SRCS:.cc=.o)
TEST_TARGETS = idec-utils-test 

all: $(TARGETS)

idec: src/main.o src/iDEC.o 
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFALAGS)

idec-steal: src/main.cc src/iDEC.cpp 
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFALAGS) -DUSE_SEARCH_HISTORY=1

kdtree_search_history_example: src/kdtree_search_history_example.o src/nanoflann.hpp
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFALAGS)

test: $(TEST_TARGETS)
	-./idec-utils-test

idec-utils-test: test/IdecUtilsTest.o
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFALAGS) $(LDFALAGS_TEST)

clean:
	-$(RM) $(TEST_OBJS) $(TEST_TARGETS) $(TARGETS) $(OBJS) $(TAR_OBJS)