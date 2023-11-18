CC=gcc
CFLAGS=-std=c99 -O3
 

ci: cInsertion.c coordReader.c
	$(CC) $(CFLAGS) cInsertion.c coordReader.c -o ci.exe -lm -fopenmp

fi: fInsertion.c coordReader.c	
	$(CC) $(CFLAGS) fInsertion.c coordReader.c -o fi.exe -lm -fopenmp

comp: ompcInsertion.c coordReader.c	
	$(CC) $(CFLAGS) ompcInsertion.c coordReader.c -o comp.exe -lm -fopenmp

fomp: ompfInsertion.c coordReader.c	
	$(CC) $(CFLAGS) ompfInsertion.c coordReader.c -o fomp.exe -lm -fopenmp 

icomp: ompcInsertion.c coordReader.c	
	icc $(CFLAGS) ompcInsertion.c coordReader.c -o icomp.exe -lm -fopenmp

ifomp: ompfInsertion.c coordReader.c	
	icc $(CFLAGS) ompfInsertion.c coordReader.c -o ifomp.exe -lm -fopenmp



