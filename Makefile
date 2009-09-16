#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE	:= sequest
KERNEL		:= kernels

HSMAIN		:= src/Main.hs
CHSFILES	:= src/Kernels.chs

SUBDIRS		:= cuda
EXTRALIBS	 = stdc++ $(KERNEL)$(LIBSUFFIX)


# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
include common.mk

