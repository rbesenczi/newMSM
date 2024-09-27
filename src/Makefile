include ${FSLCONFDIR}/default.mk

PROJNAME = newMSM
XFILES   = newmsm
LIBS     = -lfsl-newmeshreg -lfsl-newresampler -lfsl-utils
USRCXXFLAGS = -fopenmp
USRINCFLAGS += -I../../msm-newmeshreg/include/

all: ${XFILES}

%: %.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
