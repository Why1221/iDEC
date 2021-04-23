CXX = g++
# disable all warnings, as there are too many 
# please remove -DDEBUG=1 when using in practice
# CXXFLAGS = -std=c++11 -O3 -w -DDEBUG=1
CXXFLAGS = -std=c++11 -O3 -w 
LDFLAGS =

TARGETS = srs srs-compact

INT_SRC = $(wildcard *.cpp)
COMMON_INT_OBJ = $(filter-out srs.o srs-compact.o, $(INT_SRC:.cpp=.o))

TARGET_SRC = $(wildcard *.cc)
TAR_OBJ = $(TARGET_SRC:.cc=.o)


all: $(TARGETS)


srs: main.o srs.o $(COMMON_INT_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

srs-compact: main-compact.o srs-compact.o $(COMMON_INT_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

.PHONY: clean
clean:
	rm -rf $(COMMON_INT_OBJ) $(TAR_OBJ) $(TARGETS) srs.o srs-compact.o