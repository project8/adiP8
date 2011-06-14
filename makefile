CXXFLAGS := $(CXXFLAGS) -I/usr/include/ -I$(ROOTSYS)/include/ -I./ -g

TARGETDIR := .
OBJDIR := .

OBJECTS := radiation.o frequency.o paramanage.o mag_pa_tool.o fft_fitter.o bfield_ferenc.o matrix_tool.o vector_tool.o math_tool.o sim_core.o sim_pilot.o el_pa_tool.o sim_help.o sim_scatter.o eH2.o
TARGETS := adipark adi2fft adifilter adiplot

vpath %.o $(OBJDIR)


LIBS :=
LINKEDLIBS := $(shell $(ROOTSYS)/bin/root-config --noauxlibs --glibs) -lTreePlayer -lMinuit -lRootAuth -lSpectrum -lSpectrumPainter -L/usr/local/lib -lfftw3 -L/usr/lib -lm
#LINKEDLIBS := -lfftw3f

OBJECTS_WDIR = $(patsubst %.o,$(OBJDIR)/%.o,$(OBJECTS))
TARGETS_WDIR = $(patsubst %,$(TARGETDIR)/%,$(TARGETS))
TARGETOBJECTS_WDIR = $(patsubst %,$(OBJDIR)/%.o,$(TARGETS))
LIBS_WDIR = $(patsubst %,$(OBJDIR)/%,$(LIBS))

all: $(OBJECTS_WDIR) $(TARGETS_WDIR)

$(OBJECTS_WDIR) : $(OBJDIR)/%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETOBJECTS_WDIR) : $(OBJDIR)/%.o : %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGETS_WDIR) : $(TARGETDIR)/% : $(OBJDIR)/%.o $(OBJECTS_WDIR) $(LIBS_WDIR)
	$(CXX) $(CXXFLAGS) $< $(OBJECTS_WDIR) $(LIBS_WDIR) $(LINKEDLIBS) -o $@

clean:
	rm -f $(OBJECTS_WDIR)
	rm -f $(TARGETOBJECTS_WDIR)
	rm -f $(TARGETS_WDIR)
