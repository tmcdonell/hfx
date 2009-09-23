#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
#STATIC_LIB	:= libkernels.a
DYNAMIC_LIB	:= libkernels.dylib

CU_DEPS		:= $(wildcard cudpp/*) \
		   kernels.h mass.h operator.h shared_mem.h utils.h
CUFILES		:= \
		   fold.cu \
		   map.cu \
		   permute.cu \
		   radix_sort.cu \
		   replicate.cu \
		   scan.cu \
		   segmented_scan.cu \
		   zipWith.cu

CUFILES_sm_11	:= ion_series.cu

C_DEPS		:= utils.h
CFILES		:=
CCFILES		:=

INCLUDES	:=

USECUDPP	:= 1


# ------------------------------------------------------------------------------
# Rules and targets
# ------------------------------------------------------------------------------
SRCDIR		:= .
DISTROOT	:= ../dist

include ../common.mk

