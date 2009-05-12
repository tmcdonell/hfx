#
# A recipe for baking!
#

BUILD_ROOT	 = .
BUILD_TYPE	 = $(shell uname -ps | tr "[:upper:], " "[:lower:],-")
SRC_ROOT	 = $(BUILD_ROOT)/src
OUTDIR		 = $(BUILD_ROOT)/build
MAKEFLAGS	+= -I$(SRC_ROOT) --no-print-directory

TARGET		 = sequest
TARGET_P	 = sequest.p
SRC_MAIN	 = $(SRC_ROOT)/Main.hs

GHC		 = $(shell which ghc)
GHC_PACKAGES	 = 
GHC_FLAGS	 = --make -Wall -fglasgow-exts -O2
GHC_FLAGS_PROF	 = -prof -auto-all -fhpc

GHC_FLAGS	+= $(GHC_PACKAGES:%=-package %) -hidir $(OUTDIR) -odir $(OUTDIR) -i$(SRC_ROOT)
GHC_FLAGS_PROF  += -osuf p.o -hisuf p.hi

# ------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------
.PHONY: all clean
all:
	@[ -d $(OUTDIR) ] || mkdir -p $(OUTDIR)
	$(GHC) $(GHC_FLAGS) $(SRC_MAIN) -o $(TARGET)
	$(GHC) $(GHC_FLAGS) $(GHC_FLAGS_PROF) $(SRC_MAIN) -o $(TARGET_P)

clean:
	$(RM) -r $(OUTDIR) $(TARGET) $(TARGET_P)
	$(RM) -r .hpc $(TARGET_P).{aux,hp,prof,ps,tix}

