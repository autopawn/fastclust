compile:
	mkdir -p bin
	gcc -g -O4 -Wall src/*.c -lm -o bin/experiment
