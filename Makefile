FLAGS = -Wall -g -O4 -Wall

.PHONY: compile

compile:
	mkdir -p bin
	gcc $(FLAGS) sdbs/sdbs.c experiments/vectors.c -I sdbs -lm -o bin/vectors
	gcc $(FLAGS) sdbs/sdbs.c experiments/vectorsets.c -I sdbs -lm -o bin/vectorsets
