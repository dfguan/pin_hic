CC      =  gcc
#CFLAGS  =  -g -Wall -D VERBOSE -D PRINT_COVERAGE #-O2  
CFLAGS  =  -O2 -Wall -D VERBOSE -D PRINT_COVERAGE #-O2  
LDFLAGS = -lz -lm

#OBJS = gfa.o opt.o paf.o sdict.o eg.o 
PROG = pin_hic  pin_hic_it # pin_ld pin_ld_it # stat_hic #pin_hic2

.SUFFIXS:.c .o

all:$(PROG)
	mkdir -p ../bin && cp $(PROG) ../bin

pin_hic_it: bamlite.o bed.o cdict.o graph.o pin_hic_it.o sdict.o col_hic_lnks.o build_graph.o get_seq.o break_pins.o asset.o mst.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

pin_hic: bamlite.o bed.o cdict.o graph.o pin_hic.o sdict.o col_hic_lnks.o build_graph.o break_pins.o get_seq.o asset.o mst.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(PROG) ../bin

	


