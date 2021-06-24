all:	main
main:	main.o weightedGraphlib.o
		gcc -g -o main main.o weightedGraphlib.o libfdr.a -lm
main.o:	main.c 
	gcc -g -I/ -c main.c -lm
weightedGraphlib.o:	weightedGraphlib.c
	gcc -g -I/ -c weightedGraphlib.c -lm
clean:	
	rm -f *.o *~
