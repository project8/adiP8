# $Id: Makefile 230 2007-07-26 11:41:29Z s_voec01 $

OPT = -O2
CFLAGS = -W -Wall $(OPT)
LDFLAGS = -lm 

PREFIX = /usr/local
LIBDIR = $(PREFIX)/lib
BINDIR = $(PREFIX)/bin
INCLUDEDIR = $(PREFIX)/include

VERSION = 1.1.0
CURRENT = 0
CC = gcc

.PHONY: all clean dist library install

all: magmain
library: magmain_dynamic

install: library
	libtool --mode=install install -s libmagfield2.la $(LIBDIR)
	libtool --mode=install install -s magmain_dynamic $(BINDIR)/magmain
	install -m 644 magfield2.h $(INCLUDEDIR)

clean:
	rm -f *.o *.lo libmagfield2.la
	rm -f magmain magmain_dynamic
	rm -rf .libs

dist:
	mkdir magfield2-$(VERSION)
	cp Makefile *.c *.h magfield2.pdf README NEWS magfield2-$(VERSION)/
	tar cjf magfield2-$(VERSION).tar.bz2 magfield2-$(VERSION)/
	rm -rf magfield2-$(VERSION)/

libmagfield2.la: magfield2.lo
	libtool --mode=link $(CC) $(LDFLAGS) $^ -o $@ -rpath $(LIBDIR) -version-info $(CURRENT)
magmain_dynamic: magmain.o libmagfield2.la
	libtool --mode=link $(CC) $(LDFLAGS) $^ -o $@

magmain: magmain.o magfield2.o

magfield2.lo: magfield2.c magfield2.h
magfield2.o: magfield2.c magfield2.h
magmain.o: magmain.c magfield2.h

%.lo: %.c
	libtool --mode=compile $(CC) $(CFLAGS) -c $<

