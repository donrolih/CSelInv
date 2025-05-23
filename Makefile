include Makefile.in

# Define installation directory
PREFIX ?= $(HOME)/opt/CSelInv

.PHONY: all clean install

all: prepare
	$(MAKE) -C build/LIB lib
	$(MAKE) -C build/EXAMPLES
	$(MAKE) -C build/UTILITIES

prepare:
	cp -r LIB build/LIB
	cp -r EXAMPLES build/EXAMPLES
	cp -r UTILITIES build/UTILITIES
	cp Makefile.in build/

clean:
	rm -rf build/*

install:
	mkdir -p $(PREFIX)/lib
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/include
	cp build/LIB/libcsupldlt.a $(PREFIX)/lib/
	cp build/EXAMPLES/*.x $(PREFIX)/bin/
	cp build/UTILITIES/*.x $(PREFIX)/bin/