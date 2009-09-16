#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE	:= sequest
KERNEL		:= sequest-kernels

HSMAIN		:= src/Main.hs
CHSFILES	:= src/IonSeries.chs \
		   src/Sequest.chs

SUBDIRS		:= kernels
EXTRALIBS	 = stdc++ $(KERNEL)$(LIBSUFFIX)

USECUBLAS	:= 1


# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
include common.mk

