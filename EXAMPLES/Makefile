include ../Makefile.in

SUPLDLT = ../LIB/libcsupldlt.a 

LIBS = $(METIS) $(SUPLDLT) $(LAPACK) $(BLAS) -lm 

OBJS = C2Finterface.o nd2d.o 

all: clap2d5pt.x cselinv.x

.c.o:
	$(CC) $(CFLAGS)  -c $<

.f.o:
	${FC} ${FFLAGS} -c $<

cselinv.x: cselinv.o $(OBJS) $(SUPLDLT) $(METIS)
	$(LD) $(LDFLAGS) cselinv.o $(OBJS) $(LIBS) -o $@

clap2d5pt.x: clap2d5pt.o $(OBJS) $(SUPLDLT) $(METIS)
	$(LD) $(LDFLAGS) clap2d5pt.o $(OBJS) $(LIBS) -o $@


$(SUPLDLT):
	(cd ../LIB; make lib)

clean:
	rm -f clap2d5pt.x cselinv.x *.o
