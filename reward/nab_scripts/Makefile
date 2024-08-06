include config.h

SUFFIXES:
SUFFIXES: .c .cu .o

run_xmin_nmode:  xmin nmode

# Ignore for now
run_xmin:  xmin
	./xmin in.prmtop in.rst7 out.rst7 > xmin.out

run_nmode: nmode
	./nmode in.prmtop in.rst7 > nmode.out

LIBS = -L$(LIBDIR) $(FLIBS) $(NETCDFLIBF)  $(LM) 

#---------------------------------------------------------------------------
xmin:   xmin.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $< $(LIBS) $(FLIBS_ARCH)

nmode:  nmode.o
	@echo "[NABC] CC $@"
	$(VB)$(CC) -o $@ $<  $(LIBS) $(FLIBS_ARCH)

#---------------------------------------------------------------------------
clean:
	/bin/rm -f *.o xmin nmode vecs

clean_all:
	/bin/rm -f *.o xmin nmode vecs xmin.out xmode.out
