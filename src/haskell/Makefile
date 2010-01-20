#
# Baking!
#

# ------------------------------------------------------------------------------
# Input files
# ------------------------------------------------------------------------------
EXECUTABLE      := sequest
HSMAIN          := Main.hs

USECUDPP        := 1
ALGORITHMS      := smvm sort
EXTRALIBS        = stdc++ $(addsuffix $(LIBSUFFIX),$(ALGORITHMS))

# ------------------------------------------------------------------------------
# Haskell/CUDA build system
# ------------------------------------------------------------------------------
ROOTDIR := .
include ../common/common.mk
