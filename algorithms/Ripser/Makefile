build: ripser


all: ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug


ripser: main.cpp ripser.cpp
	c++ -std=c++11 main.cpp ripser.cpp -o ripser -Ofast -D NDEBUG

ripser-coeff: main.cpp ripser.cpp
	c++ -std=c++11 main.cpp ripser.cpp -o ripser-coeff -Ofast -D NDEBUG -D USE_COEFFICIENTS

ripser-reduction: main.cpp ripser.cpp
	c++ -std=c++11 main.cpp ripser.cpp -o ripser-reduction -Ofast -D NDEBUG -D ASSEMBLE_REDUCTION_MATRIX

ripser-coeff-reduction: main.cpp ripser.cpp
	c++ -std=c++11 main.cpp ripser.cpp -o ripser-coeff-reduction -Ofast -D NDEBUG -D USE_COEFFICIENTS -D ASSEMBLE_REDUCTION_MATRIX

ripser-debug: main.cpp ripser.cpp
	c++ -std=c++11 main.cpp ripser.cpp -o ripser-debug -g


clean:
	rm -f ripser ripser-coeff ripser-reduction ripser-coeff-reduction ripser-debug
