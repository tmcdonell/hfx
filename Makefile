#
#
#

BUILD_ROOT	 = .
BUILD_TYPE	 = $(shell uname -ps | tr "[:upper:], " "[:lower:],-")
SRC_ROOT	 = $(BUILD_ROOT)/src
OUTDIR		 = $(BUILD_ROOT)/build
MAKEFLAGS	+= -I$(SRC_ROOT) --no-print-directory

TARGET		 = sequest
SRC_MAIN	 = $(SRC_ROOT)/Main.hs

GHC		 = $(shell which ghc)
GHC_PACKAGES	 = 
GHC_FLAGS	 = -Wall --make -fglasgow-exts -O2
GHC_FLAGS	+= -prof -auto-all -fhpc

# ------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------
.PHONY: all clean
all:
	@[ -d $(OUTDIR) ] || mkdir -p $(OUTDIR)
	$(GHC) $(GHC_FLAGS) $(GHC_PACKAGES:%=-package %) -hidir $(OUTDIR) -odir $(OUTDIR) -i$(SRC_ROOT) $(SRC_MAIN) -o $(TARGET)

clean:
	$(RM) -r $(OUTDIR) $(TARGET)
	$(RM) -r .hpc $(TARGET).{aux,hp,prof,ps,tix}

