compile:
	mkdir -p bin
	gcc -Wall -g -O4 -Wall src/*.c -lm -o bin/experiment
