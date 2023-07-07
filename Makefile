CC=clang++
CFLAGS= -Wall -Wextra -O3 -std=c++20


lessgex: lessgex.cpp lessgex.h
	$(CC) $(CFLAGS) -o lessgex lessgex.cpp

clean:
	rm -f lessgex
	rm -f lessgex.o