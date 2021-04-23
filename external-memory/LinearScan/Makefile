CXX = g++
CXXFLAGS = -std=c++14 -O3 -Wall -I/usr/include/hdf5/serial -I../Common
LDFLAGS = -L /usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5
RM = gio trash -f


TARGETS = linear-scan-compact-example   linear-scan linear-scan-l1

INT_SRC = $(wildcard *.cpp)
INT_OBJ = $(INT_SRC:.cpp=.o)

TARGET_SRC = $(wildcard *.cc)
TAR_OBJ = $(TARGET_SRC:.cc=.o)


all: $(TARGETS)


linear-scan-compact-example: linear_scan_compact_example.o $(INT_OBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

linear-scan: linear_scan_main.o $(INT_OBJ) Hdf5File.h
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

linear-scan-l1: linear_scan_main_l1.o $(INT_OBJ) Hdf5File.h
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)


.PHONY: clean
clean:
	$(RM) $(INT_OBJ) $(TAR_OBJ) $(TARGETS)

distclean: clean
	$(RM) linear-scan.txt demo-results.txt gnd.txt
