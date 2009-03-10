#
#
#

BUILD_ROOT	 = .
SRC_ROOT	 = $(BUILD_ROOT)/src
OUTDIR		 = $(BUILD_ROOT)/build-$(shell uname -ps | tr "[:upper:], " "[:lower:],-")
MAKEFLAGS	+= -I$(SRC_ROOT) --no-print-directory

TARGET		 = sequest
SRC_MAIN	 = $(SRC_ROOT)/sequest.hs

GHC		 = $(shell which ghc)
GHC_PACKAGES	 =
GHC_FLAGS	 = --make -Wall -fglasgow-exts -hidir $(OUTDIR) -odir $(OUTDIR)

# ------------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------------
all:
	@[ -d $(OUTDIR) ] || mkdir -p $(OUTDIR)
	$(GHC) $(GHC_FLAGS) -i$(SRC_ROOT) $(SRC_MAIN) -o $(TARGET)

.PHONY: clean
clean:
	$(RM) -r $(OUTDIR) $(TARGET)

