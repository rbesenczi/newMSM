include ${FSLCONFDIR}/default.mk

PROJNAME = newresampler-demo
XFILES   = metric-resample surface-resample smoothing NN-resample
LIBS     = -lfsl-newresampler -lboost_program_options -lfsl-miscmaths
USRCXXFLAGS += -Wall -fopenmp

all: ${XFILES}

%: %.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
