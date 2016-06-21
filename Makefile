CC=g++
#CFLAGS=-Wall -Werror -shared -fPIC 
#CFLAGS=-Wall -shared -fPIC 
#$(CC) -shared -Wl,soname,libputN.so -o libputN.so libputN.o
	#$(CC) --whole-archive -shared -Wl,-soname,libputN.so -o libputN.so libputN.o
	#$(CC) -O2 -Wall -fPIC -lstdc++ -g -I. -o libputN.o putN.cpp
	#ld -G -o libputN.so putN.o

all:  putN.cpp 
	$(CC) -O2 -fPIC -c putN.cpp
	$(CC) -shared -Wl,-soname,libputN.so putN.o -lsupc++ -o libputN.so

clean:
	rm libputN.so
	rm putN.o
