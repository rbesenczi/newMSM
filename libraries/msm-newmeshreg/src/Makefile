include ${FSLCONFDIR}/default.mk

PROJNAME    = NewMeshReg
SOFILES     = libfsl-newmeshreg.so
LIBS        = -lfsl-newresampler -lfsl-newimage -lfsl-miscmaths -lfsl-utils
USRCXXFLAGS = -fopenmp
USRINCFLAGS += -DHAS_FPD -DHAS_HOCR -I../include/

OBJS       = reg_tools.o featurespace.o similarities.o DiscreteCostFunction.o \
             DiscreteModel.o mesh_registration.o rigid_costfunction.o \
             group_mesh_registration.o DiscreteGroupModel.o DiscreteGroupCostFunction.o

all: ${SOFILES}

libfsl-newmeshreg.so: ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}
