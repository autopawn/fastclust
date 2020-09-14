FLAGS = -Wall -g -O4 -Wall

compile:
	mkdir -p bin
	gcc $(FLAGS) sdbs/*.c -lm -o bin/experiment
