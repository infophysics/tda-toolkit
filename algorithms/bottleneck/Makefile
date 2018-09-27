CC=g++
INC=/Applications/MATLAB_R2016b.app/extern/include
CFLAGS=-Wall -O3 -I$(INC)
EXECUTABLES=bottleneck bottleneckC

all:
	$(CC) $(CFLAGS) bottleneckDist.cpp -o bottleneck
clean:
	rm -rf $(EXECUTABLES)
