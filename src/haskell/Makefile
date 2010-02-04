#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE      := sequest
HSMAIN          := Main.hs

ALGORITHMS      := enum filter ion_series mvm permute radix_sort replicate scan smvm
EXTRALIBS       += stdc++ $(addsuffix $(LIBSUFFIX),$(ALGORITHMS))

# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
ROOTDIR := .
include ../common/common.mk
