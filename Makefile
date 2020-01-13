
rsaJR:rsaJR.o main.o
	gcc -o rsaJR rsaJR.o main.o

rsaJR_crible:rsaJR_crible.o main.o
	gcc -o rsaJR_crible rsaJR_crible.o main.o

rsaJR_crible.o: rsaJR_crible.c rsaJR.h
	gcc -o rsaJR_crible.o rsaJR_crible.c -W -Wall -ansi -pedantic

rsaJR.o:rsaJR.c rsaJR.h
	gcc -o rsaJR.o rsaJR.c -W -Wall -ansi -pedantic

main.o:main.c 
	gcc -o main.o -c main.c -W -Wall -ansi -pedantic

clean:
	rm -rf *.o

mrproper: clean
	rm -rf hello
