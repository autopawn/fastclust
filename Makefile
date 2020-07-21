compile:
	mkdir -p bin
	gcc -Wall -g -O4 -D DEBUGINFO -Wall src/*.c -lm -o bin/experiment
