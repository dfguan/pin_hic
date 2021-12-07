CC      =  gcc
CFLAGS  =  -g -Wall -D VERBOSE -D PRINT_COVERAGE -I ../src#-O2  
#CFLAGS  =  -O2 -Wall -D VERBOSE -D PRINT_COVERAGE -I ../src #-O2  
LDFLAGS = -lz -lm

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
UTIL = bam2hig_bed 

.SUFFIXS:.c .o

all:$(UTIL)
	

bam2hig_bed:bam2hig_bed.o ../src/bamlite.o ../src/cdict.o ../src/graph.o   ../src/asset.o ../src/sdict.o ../src/bed.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(UTIL) 

	


