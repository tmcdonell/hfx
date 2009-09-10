#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
#STATIC_LIB	:= libsequest-kernels.a
DYNAMIC_LIB	:= libsequest-kernels.dylib

CU_DEPS		:= $(wildcard cudpp/*) \
		   kernels.h mass.h operator.h shared_mem.h utils.h \
		   map.cuh reduce.cuh scan.cuh zipWith.cuh
CUFILES		:=
CUFILES_sm_11	:= prelude.cu ion_series.cu

C_DEPS		:= utils.h
CFILES		:=
CCFILES		:=

INCLUDES	:=


# ------------------------------------------------------------------------------
# Rules and targets
# ------------------------------------------------------------------------------
SRCDIR		:= .
DISTROOT	:= ../dist

include ../common.mk

