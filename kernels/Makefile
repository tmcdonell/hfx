#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
#STATIC_LIB	:= libsequest-kernels.a
DYNAMIC_LIB	:= libsequest-kernels.dylib

CU_DEPS		:= $(wildcard cudpp/*) \
		   kernels.h mass.h operator.h fold.cuh scan.cuh zipWith.cuh
CUFILES		:=
CUFILES_sm_11	:= test.cu

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

