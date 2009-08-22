#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE	:= sequest
KERNEL		:= sequest-kernels

HSMAIN		:= src/Main.hs
SUBDIRS		:= kernels
EXTRALIBS	:= $(KERNEL)


# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
include common.mk

