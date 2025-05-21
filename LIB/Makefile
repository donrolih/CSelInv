include ../Makefile.in

# LLIN: Add -D__DEBUG for verbose information

.F.o:
	$(FC) $(FFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

OBJS =assmb.o \
btree2.o \
betree.o \
bfinit.o \
blkfct.o \
blkfc2.o \
blkslv.o \
chordr.o \
chlsup.o \
create.o \
epost2.o \
etordr.o \
etpost.o \
etree.o \
fcnthn.o \
fsup1.o \
fsup2.o \
fnsplt.o \
fntsiz.o \
genmmd.o \
getrhs.o \
gtimer.o \
igathr.o \
inpnv.o \
invinv.o \
ldindx.o \
lstats.o \
mmdint.o \
mmdelm.o \
mmdupd.o \
mmdnum.o \
mmpy.o \
mmpy1.o \
mmpy2.o \
mmpy4.o \
mmpy8.o \
mmpyi.o \
ordmmd.o \
ordnat.o \
pchol.o \
sfinit.o \
smxpy1.o \
smxpy2.o \
smxpy4.o \
smxpy8.o \
symfct.o \
symfc2.o \
flo2ho.o \
ilo2ho.o \
dreadhb.o \
dsymmv.o \
exdiag.o \
exdiagblk.o \
selinvblk.o \
sortix.o \
dumpL.o \
readcmatrix.o \
readdmatrix.o \
readmatrixheader.o\
zdotu.o\
dscal1.o 

# Comment:
# dscale1 is renamed to avoid conflict with BLAS
# zdotu is needed on some machinetypes such as OSX

lib: $(OBJS)
	$(AR) libcsupldlt.a $(OBJS)
	$(RANLIB) libcsupldlt.a

clean:
	rm -f *.o *.a
