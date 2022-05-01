CC =g++
TARGET =main-test

SOURCES =GVCFWriter.cpp main_test.cpp

CPPFLAGS=-std=c++17 -Wall -Wextra -g -O0 -march=native -I ./htslib -L ./
LDFLAGS=
LDLIBS=-L. libhts.so.3

# define list of objects
# OBJSC=$(SOURCES:.cpp=.o)
OBJS=GVCFWriter.o main_test.o

# the target is obtained linking all .o files
all: $(SOURCES) $(TARGET)

$(TARGET): $(OBJS) sample.o
	$(CC) $(OBJS) sample.o -o $(TARGET) $(LDFLAGS) $(LDLIBS)

purge: clean

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)