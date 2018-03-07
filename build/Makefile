# makefile
#
#
#  Directives :
#
#    

SHELL=/bin/bash

include envfile
export $(shell sed 's/=.*//' envfile)

ifeq (${FC},gfortran)

  OPT=-O3
  ifeq (${OMP},yes)
    OPT+=-fopenmp
  endif

endif

# determine the main directory for this experiment
GeoBarSDIR=$(shell pwd | sed 's/\GeoBarS.*/GeoBarS/')
SRCDIR=${GeoBarSDIR}/src/


# ---------------------------------------------------------------------- #
#                           ~/src/common/                                #
# ---------------------------------------------------------------------- #

# Modules
ModelPrecision.o :
	${FC} ${OPT} -c ${SRCDIR}common/ModelPrecision.f90 -o $@

ConstantsDictionary.o : ModelPrecision.o
	${FC} ${OPT} -c ${SRCDIR}common/ConstantsDictionary.f90 -o $@

LinkedList_Class.o :
	${FC} ${OPT} -c ${SRCDIR}common/LinkedList_Class.f90 -o $@

KeyRing_Class.o : CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/KeyRing_Class.f90 -o $@

HashTable_Class.o : ConstantsDictionary.o LinkedList_Class.o
	${FC} ${OPT} -c ${SRCDIR}common/HashTable_Class.f90 -o $@

CommonRoutines.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}common/CommonRoutines.f90 -o $@

Timing.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}common/Timing.f90 -o $@

# Test programs
   
# ---------------------------------------------------------------------- #
#                           ~/src/spectralops/                           #
# ---------------------------------------------------------------------- #

# Modules
Quadrature.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/Quadrature.f90 -o $@

InterpolationSupportRoutines.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o 
	${FC} ${OPT} -c ${SRCDIR}spectralops/InterpolationSupportRoutines.f90 -o $@

Lagrange_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o InterpolationSupportRoutines.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/Lagrange_Class.f90 -o $@

NodalStorage_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                       Quadrature.o Lagrange_Class.o
	${FC} ${OPT} -c ${SRCDIR}spectralops/NodalStorage_Class.f90 -o $@



# ---------------------------------------------------------------------- #
#                            ~/src/geom/                                 #
# ---------------------------------------------------------------------- #

GOBJS = ModelPrecision.o \
        ConstantsDictionary.o \
        CommonRoutines.o \
        LinkedList_Class.o \
        HashTable_Class.o \
        Timing.o \
        Quadrature.o \
        InterpolationSupportRoutines.o \
        Lagrange_Class.o \
        NodalStorage_Class.o \
        DGSEMSolutionStorage_2D_Class.o \
        Node_Class.o \
        Edge_Class.o \
        Curve_Class.o \
        MappedGeometry_2D_Class.o \
        QuadElement_Class.o \
        QuadMesh_Class.o \
        Params_Class.o \
        BoundaryCommunicator_Class.o

         
GeneratePCMesh_2D : ${GOBJS} GenerateMeshFiles_2D.o 
	${FC} ${OPT} ${GOBJS} GenerateMeshFiles_2D.o  -o $@


## Programs
GenerateMeshFiles_2D.o : ModelPrecision.o NodalStorage_Class.o QuadMesh_Class.o \
                         BoundaryCommunicator_Class.o Params_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/GenerateMeshFiles_2D.f90 -o $@

## Modules 

GeometryParams_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o
	${FC} ${OPT} -c ${SRCDIR}geom/GeometryParams_Class.f90 -o $@

Node_Class.o : ModelPrecision.o  ConstantsDictionary.o  LinkedList_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/Node_Class.f90 -o $@

Edge_Class.o : ModelPrecision.o  ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/Edge_Class.f90 -o $@

QuadElement_Class.o : ModelPrecision.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadElement_Class.f90 -o $@

Curve_Class.o : ModelPrecision.o ConstantsDictionary.o Lagrange_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/Curve_Class.f90 -o $@

MappedGeometry_2D_Class.o : ModelPrecision.o ConstantsDictionary.o  CommonRoutines.o \
                            Lagrange_Class.o Curve_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/MappedGeometry_2D_Class.f90 -o $@

QuadMesh_Class.o : ModelPrecision.o ConstantsDictionary.o  LinkedList_Class.o HashTable_Class.o \
                   Quadrature.o Curve_Class.o MappedGeometry_2D_Class.o QuadElement_Class.o \
                   Edge_Class.o Node_Class.o
	${FC} ${OPT} -c ${SRCDIR}geom/QuadMesh_Class.f90 -o $@

## ---------------------------------------------------------------------- #
##                            ~/src/highend/                              #
## ---------------------------------------------------------------------- #


ifdef MODDIR
GeoBarSParams_Class.o : ModelPrecision.o CommonRoutines.o ConstantsDictionary.o
	${FC} ${OPT} -c ${MODDIR}/GeoBarSParams_Class.f90 -o $@

GeoBarS_Driver.o : ModelPrecision.o ConstantsDictionary.o QuadMesh_Class.o GeoBarS_Class.o
	${FC} ${OPT} -c ${MODDIR}/GeoBarS_Driver.f90 -o $@
else
GeoBarSParams_Class.o : ModelPrecision.o CommonRoutines.o ConstantsDictionary.o
	${FC} ${OPT} -c ${SRCDIR}highend/GeoBarSParams_Class.f90 -o $@

GeoBarS_Driver.o : ModelPrecision.o ConstantsDictionary.o QuadMesh_Class.o GeoBarS_Class.o
	${FC} ${OPT} -c ${SRCDIR}highend/GeoBarS_Driver.f90 -o $@
endif

GeoBarS_Class.o : ModelPrecision.o ConstantsDictionary.o CommonRoutines.o NodalStorage_Class.o \
                  QuadMesh_Class.o Edge_Class.o Node_Class.o MappedGeometry_2D_Class.o GeoBarSParams_Class.o
	${FC} ${OPT} -c ${SRCDIR}highend/GeoBarS_Class.f90 -o $@

GEOBJ= ModelPrecision.o \
       ConstantsDictionary.o \
       LinkedList_Class.o \
       HashTable_Class.o \
       Quadrature.o \
       CommonRoutines.o \
       InterpolationSupportRoutines.o \
       Lagrange_Class.o \
       Curve_Class.o \
       MappedGeometry_2D_Class.o \
       QuadElement_Class.o \
       Edge_Class.o \
       Node_Class.o \
       QuadMesh_Class.o \
       NodalStorage_Class.o \
       GeoBarSParams_Class.o \
       GeoBarS_Class.o \
       GeoBarS_Driver.o

geobars : ${GEOBJ}
	${FC} ${OPT} ${GEOBJ} -o $@

# ---------------------------------------------------------------------- #
#                          Other directives                              #
# ---------------------------------------------------------------------- #
.PHONY : clean

clean :
	rm *.o *.mod



















