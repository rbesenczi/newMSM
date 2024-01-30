include ${FSLCONFDIR}/default.mk

PROJNAME = newresampler-demo
XFILES   = metric-resample surface-resample smoothing NN-resample applywarp
LIBS     = -lfsl-newresampler -lboost_program_options -lfsl-miscmaths
USRCXXFLAGS += -fopenmp

all: ${XFILES}

%: %.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
