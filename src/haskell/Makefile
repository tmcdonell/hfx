#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE      := sequest
HSMAIN          := Main.hs

ALGORITHMS      := filter ion_series mvm permute radix_sort scan
EXTRALIBS       += stdc++ $(addsuffix $(LIBSUFFIX),$(ALGORITHMS))

# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
ROOTDIR := .
include ../common/common.mk
