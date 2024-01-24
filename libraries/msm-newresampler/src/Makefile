include ${FSLCONFDIR}/default.mk

PROJNAME = newresampler
SOFILES  = libfsl-newresampler.so
LIBS     = -lfsl-newmeshclass -lfsl-miscmaths
OBJS     = meshException.o point.o mpoint.o triangle.o\
		   mesh.o node.o octree.o resampler.o
USRCXXFLAGS += -fopenmp

all: libfsl-newresampler.so

libfsl-newresampler.so: ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}
