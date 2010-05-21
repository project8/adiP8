.SUFFIXES: .c . .o

HEADERS = -I$(ROOTSYS)/include/  -I/usr/local/include 
LIBS = -L$(ROOTSYS)/lib -lCore -lCint -lMatrix -lTreePlayer -lMinuit -lHist -lNet -lGraf -lTree -lGpad -lPostscript -lPhysics -lSpectrum -lMathCore -lSpectrumPainter -lRIO -L/usr/local/lib -lfftw3 -lgsl -L/usr/lib -lm


GCC_OPT=-O1 -funroll-loops -ffast-math -Wall 
# LINK_OPT=-L./local/lib/ -I./ 
# modified by bmonreal
LINK_OPT= -I./

.c: 
	g++ $(GCC_OPT) $(LINK_OPT) -s -o $* $*.c -lm -lcurses
.c.o: 
	g++  $(GCC_OPT) -c -o $*.o $*.c 


fft: adi2fft.o fft_fitter.o radiation.o paramanage.o mag_pa_tool.o frequency.o bfield_ferenc.o matrix_tool.o vector_tool.o math_tool.o
	g++  $(LINK_OPT) -o adi2fft adi2fft.o fft_fitter.o radiation.o paramanage.o mag_pa_tool.o frequency.o bfield_ferenc.o matrix_tool.o vector_tool.o math_tool.o $(LIBS)

autogen_runfile: autogen_runfile.c
	g++ -o autogen_runfile autogen_runfile.c $(HEADERS) $(LIBS)

radiation.o: radiation.c frequency.h
	g++  -c -o radiation.o  radiation.c $(HEADERS) 

adi2fft.o: adi2fft.c paramanage.h fft_fitter.h radiation.h mag_pa_tool.h frequency.h
	g++  -c -o adi2fft.o  adi2fft.c $(HEADERS)

fft_fitter.o: fft_fitter.c fft_fitter.h
	g++  -c -o fft_fitter.o fft_fitter.c $(HEADERS)	

frequency.o: frequency.c frequency.h
	g++ -c -o frequency.o frequency.c $(HEADERS)

adipark: adipark.o matrix_tool.o sim_core.o vector_tool.o math_tool.o sim_pilot.o mag_pa_tool.o bfield_ferenc.o el_pa_tool.o paramanage.o sim_help.o sim_scatter.o eH2.o
	g++ $(GCC_OPT) $(LINK_OPT) -o adipark matrix_tool.o adipark.o sim_core.o vector_tool.o math_tool.o sim_pilot.o  mag_pa_tool.o bfield_ferenc.o el_pa_tool.o paramanage.o sim_help.o sim_scatter.o eH2.o -lm 

adipark.o: adipark.c el_pa_tool.h mag_pa_tool.h paramanage.h
	g++ $(GCC_OPT) -c -o adipark.o adipark.c

sim_pilot.o: sim_pilot.c sim_pilot.h sim_core.h sim_help.h paramanage.h sim_scatter.h
	g++ $(GCC_OPT) -c -o sim_pilot.o sim_pilot.c

sim_core.o: sim_core.c sim_core.h vector_tool.h math_tool.h el_pa_tool.h mag_pa_tool.h paramanage.h sim_help.h sim_scatter.h
	g++ $(GCC_OPT) -c -o sim_core.o sim_core.c

el_pa_tool.o: el_pa_tool.c el_pa_tool.h matrix_tool.h paramanage.h
	g++ $(GCC_OPT) -c -o el_pa_tool.o el_pa_tool.c

mag_pa_tool.o: mag_pa_tool.c mag_pa_tool.h matrix_tool.h paramanage.h bfield_ferenc.h
	g++ $(GCC_OPT) -c -o mag_pa_tool.o mag_pa_tool.c

bfield_ferenc.o: bfield_ferenc.h paramanage.h
	g++ $(GCC_OPT) -c -o bfield_ferenc.o bfield_ferenc.c

sim_scatter.o: sim_scatter.c sim_scatter.h eH2.h paramanage.h
	g++ $(GCC_OPT) -c -o sim_scatter.o sim_scatter.c 

eH2.o: eH2.c eH2.h
	g++ $(GCC_OPT) -c -o eH2.o eH2.c 

paramanage.o: paramanage.h
	g++ $(GCC_OPT) -c -o paramanage.o paramanage.c

sim_help.o: sim_help.c sim_help.h
	g++ $(GCC_OPT) -c -o sim_help.o sim_help.c 

matrix_tool.o: matrix_tool.c matrix_tool.h
	g++ $(GCC_OPT) -c -o matrix_tool.o matrix_tool.c

vector_tool.o: vector_tool.c vector_tool.h math_tool.h
	g++ $(GCC_OPT) -c -o vector_tool.o vector_tool.c 

math_tool.o: math_tool.c math_tool.h
	g++ $(GCC_OPT) -c -o math_tool.o math_tool.c 


clean: clean-src
	rm -f *.o .depend $(TARGETS)

clean-src:
	rm -f *~
	rm -f *.bak
