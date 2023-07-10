CC=clang++
CFLAGS= -Wall -Wextra -O3 -std=c++20 -Wshadow -g


lessgex: lessgex.cpp lessgex.h
	$(CC) $(CFLAGS) -o lessgex lessgex.cpp

bench: bench.cpp lessgex.h
	$(CC) $(CFLAGS) -o bench bench.cpp -lpcre2-8 

clean:
	rm -f bench
	rm -f bench.o
	rm -f lessgex
	rm -f lessgex.o