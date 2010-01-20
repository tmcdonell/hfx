#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE      := sequest
HSMAIN          := Main.hs

ALGORITHMS      := radix_sort scan smvm
EXTRALIBS       += stdc++ $(addsuffix $(LIBSUFFIX),$(ALGORITHMS))

# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
ROOTDIR := .
include ../common/common.mk
