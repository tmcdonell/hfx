#
# Common bridge between Haskell and CUDA build systems
# 
# For the variables set below, it was required to edit the common.mk file from
# the SDK to be optional defines (?=).
#

# ------------------------------------------------------------------------------
# Common CUDA build system
# ------------------------------------------------------------------------------

CUDA_SDK_PATH	 = /Developer/CUDA

SRCDIR		 = src/
DISTROOT	 = dist
BINDIR		 = $(DISTROOT)/bin
ROOTOBJDIR	 = $(DISTROOT)/obj
LIBDIR		 = $(CUDA_SDK_PATH)/lib
COMMONDIR	 = $(CUDA_SDK_PATH)/common

GHC_FLAGS	 = --make -Wall -fglasgow-exts -lstdc++ \
		   -i$(SRCDIR) -i$(ROOTOBJDIR)/$(BINSUBDIR) \
		   -odir $(ROOTOBJDIR)/$(BINSUBDIR) -hidir $(ROOTOBJDIR)/$(BINSUBDIR)

# small hack: dependencies of the target, but not passed to object linker
#
CUBINS		+= $(patsubst %.chs,$(OBJDIR)/%.hs,$(notdir $(CHSFILES)))

ifeq ($(dbg),1)
    GHC_FLAGS	+= -prof -auto-all -caf-all -fhpc
else
    GHC_FLAGS	+= -O3
endif

ifeq ($(suffix $(HSMAIN)),.chs)
    LINK	 = ghc $(GHC_FLAGS) $(OBJDIR)/$(HSMAIN:%.chs=%.hs)
else
    LINK	 = ghc $(GHC_FLAGS) $(SRCDIR)$(HSMAIN)
endif

# ------------------------------------------------------------------------------
# Rules
# ------------------------------------------------------------------------------
include $(COMMONDIR)/common.mk

$(OBJDIR)/%.hs : $(SRCDIR)%.chs
	$(VERBOSE)c2hs --include=$(OBJDIR) $(INCLUDES:%=--cppopts=%) --output-dir=$(OBJDIR) --output=$(notdir $@) $<

spotless : clean
	$(VERBOSE)rm -rf $(DISTROOT)
	$(VERBOSE)rm -f  $(patsubst %.cu,%.linkinfo,$(notdir $(CUFILES)))

