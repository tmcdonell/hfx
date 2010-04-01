#
# First build all of the CUDA projects, producing a series of algorithm
# libraries. Then build and link the main project against these.
#

ALGORITHMS := $(shell find src/cuda -name Makefile)
PROJECTS   := $(ALGORITHMS)
CABAL      := cabal

ifeq ($(dbg),1)
    CABALFLAGS += -fdebug
endif
ifeq ($(verbose),1)
    CABALFLAGS += --verbose
    VERBOSE    :=
else
    VERBOSE    := @
endif


%.do :
	$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

all : $(addsuffix .do,$(PROJECTS))
	$(VERBOSE)$(CABAL) configure $(CABALFLAGS)
	$(VERBOSE)$(CABAL) build $(CABALFLAGS)
	@echo "Finished building all"

clean : $(addsuffix .do,$(PROJECTS))
	$(VERBOSE)$(CABAL) clean
	@echo "Finished cleaning all"

clobber : $(addsuffix .do,$(PROJECTS))
	$(VERBOSE)$(CABAL) clean
	@echo "Finished cleaning all"

