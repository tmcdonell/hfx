#
# First build all of the CUDA projects, producing a series of algorithm
# libraries. Then build and link the main project against these.
#

ALGORITHMS := $(shell find src/cuda -name Makefile)
PROJECTS   := $(ALGORITHMS) src/sequest/Makefile


%.do :
	$(MAKE) -C $(dir $*) $(MAKECMDGOALS)

all : $(addsuffix .do,$(PROJECTS))
	@echo "Finished building all"

clean : $(addsuffix .do,$(PROJECTS))
	@echo "Finished cleaning all"

clobber : $(addsuffix .do,$(PROJECTS))
	@echo "Finished cleaning all"

