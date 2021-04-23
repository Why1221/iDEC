CXX = g++ 
CXXFLAGS = -std=c++14 -Wall 
LDFLAGS = 
RM = gio trash -f


COMMON_HDR = Exception.h
TARGETS =  ann-result-writer-test timer-test
SRCS = $(wildcard *.cc)
OBJS = $(SRCS:.cc=.o)


all: $(TARGETS)

ann-result-writer-test: AnnResultWriterTest.o AnnResultWriter.hpp $(COMMON_HDR)
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

timer-test: TimerTest.o Timer.hpp
	-$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

.PHONY: clean 

clean:
	-$(RM) $(OBJS) $(TARGETS) AnnResultWriterTest TimerTest