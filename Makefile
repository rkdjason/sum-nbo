#Makefile

all: sum-nbo

sum-nbo: main.cpp
	g++ -o sum-nbo main.cpp
	
clean:
	rm -f sum-nbo
