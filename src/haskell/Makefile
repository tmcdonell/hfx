#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE      := sequest
HSMAIN          := Main.hs

ALGORITHMS      := filter ion_series mvm sort
EXTRALIBS       += stdc++ $(addsuffix $(LIBSUFFIX),$(ALGORITHMS))
USECUDPP	:= 0

# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
ROOTDIR := .
include ../common/common.mk
