import os

# General System:
Flags = [Split('-g -Wall -Wextra -Werror')]
IncludePath = ['/usr/include', './include', '.']
LibPath = [Split('/usr/lib')]
Libs = []

# Deal With ROOT:
rootPath = [os.environ['ROOTSYS']+'/include']
rootconf = os.popen('root-config --noauxcflags --glibs').read().strip().split()
rootLibPath = []
rootLibs = []
for rc in rootconf:
    if rc.startswith('-L'):
        rootLibPath.append(rc.replace('-L',''))
    if rc.startswith('-l'):
        rootLibs.append(rc.replace('-l',''))

# Deal with FFTW3
fftwLibs = ['fftw3','m']

# Setup build environment
#bld_rootdict = Builder(action = os.environ['ROOTSYS']+"/bin/rootcint -f $TARGET -c -g -Wall -Wextra -Werror -p $SOURCES && mv `echo $TARGET | sed -e 's/cpp/h/'`",
bld_rootdict = Builder(action = os.environ['ROOTSYS']+"/bin/rootcint -f $TARGET -c -g -Wall -Wextra -Werror -p $SOURCES && mv $TARGET `echo $TARGET | sed -e 's/include/src/'`",
                       suffix = '.cpp',
                       src_suffix = '.hpp')

env = Environment(CC='g++', CPPFLAGS = Flags, CPPPATH = IncludePath + rootPath, LIBPATH = LibPath + rootLibPath, LIBS = Libs + rootLibs + fftwLibs)
env['ENV']['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH']
env['BUILDERS']['RootDict'] = bld_rootdict

# Actually make the build rules
env.Program(target = 'bin/adipark', source = Split('src/adipark.c src/mag_pa_tool.c src/magfield3.c src/bfield_ferenc.c src/eH2.c src/math_tool.c src/vector_tool.c src/matrix_tool.c src/array.c src/sim_pilot.c src/sim_core.c src/sim_scatter.c src/sim_help.c src/paramanage.c src/el_pa_tool.c'))

env.Program(target = 'bin/adi2fft', source = Split('src/adi2fft.c src/fft_fitter.c src/paramanage.c src/radiation.c'))

env.Program(target = 'bin/adifilter', source = Split('src/adifilter.c src/paramanage.c src/radiation.c'))

env.Program(target = 'bin/adiplot', source = Split('src/adiplot.c src/paramanage.c src/radiation.c'))

env.Program(target = 'bin/magsource', source = Split('src/magsource.c src/magfield3.c src/array.c'))

env.RootDict(target = 'include/clp8.cpp', source = Split('include/parameter include/window_fft'))

env.SharedLibrary(target = 'lib/clp8', source = Split('src/clp8.cpp src/parameter.cpp src/window_fft.cpp'))#, CPPPATH = IncludePath + rootPath + ['.'])
