compile:
	mkdir -p bin
	gcc -g -Wall src/*.c -lm -o bin/experiment
