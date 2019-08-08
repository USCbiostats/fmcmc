check:
	$(MAKE) build && R CMD check ../fmcmc_*.tar.gz && $(MAKE) clean

checkv:
	$(MAKE) build && R CMD check --use-valgrind ../fmcmc_*.tar.gz && $(MAKE) clean

build:
	cd ../ && R CMD build fmcmc/ && cd fmcmc

clean:
	rm  ../fmcmc_*.tar.gz

install:
	R CMD INSTALL --preclean .

.PHONY: check checkv build clean 
