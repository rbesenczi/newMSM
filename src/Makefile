include ${FSLCONFDIR}/default.mk

PROJNAME = newMSM
XFILES   = newmsm
LIBS     = -lfsl-newmeshreg -lfsl-newresampler -lfsl-utils
USRCXXFLAGS = -fopenmp
USRINCFLAGS = -DHAS_HOCR

all: ${XFILES}

%: %.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}
