all: ga

ga: ga.cpp
	g++ -std=c++17 -o ga -O3 ga.cpp

run: ga
	./ga < maxcut.in > maxcut.out

clean:
	rm ga
