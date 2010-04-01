#
# First build all of the CUDA projects, producing a series of algorithm
# libraries. Then build and link the main project against these.
#

ALGORITHMS := $(shell find src/cuda -name Makefile)
PROJECTS   := $(ALGORITHMS)
CABAL      := cabal


%.do :
	$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

all : $(addsuffix .do,$(PROJECTS))
	@$(CABAL) configure
	@$(CABAL) build
	@echo "Finished building all"

clean : $(addsuffix .do,$(PROJECTS))
	@$(CABAL) clean
	@echo "Finished cleaning all"

clobber : $(addsuffix .do,$(PROJECTS))
	@$(CABAL) clean
	@echo "Finished cleaning all"

