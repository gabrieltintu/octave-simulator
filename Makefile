# Tintu Gabriel-Claudiu 2023 - 2024

# compiler setup
CC=gcc
CFLAGS=-Wall -Wextra -std=c99

# define targets
TARGETS=my_octave

build: $(TARGETS)

my_octave: my_octave.c functions.c
	$(CC) $(CFLAGS) my_octave.c functions.c -o my_octave

pack:
	zip -FSr octave_simulator.zip README Makefile *.c *.h

clean:
	rm -f $(TARGETS)

.PHONY: pack clean