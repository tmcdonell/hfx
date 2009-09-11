#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE	:= sequest
KERNEL		:= sequest-kernels

HSMAIN		:= src/Main.hs
CHSFILES	:= src/IonSeries.chs

SUBDIRS		:= kernels
EXTRALIBS	 = $(KERNEL)$(LIBSUFFIX)


# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
include common.mk

