
NAUTY_DIR = nauty

PYTHON = python
PIP = pip
python_version_full := $(wordlist 2,4,$(subst ., ,$(shell $(PYTHON) --version 2>&1)))
python_version_major := $(word 1,${python_version_full})
python_version_minor := $(word 2,${python_version_full})
python_version_patch := $(word 3,${python_version_full})
machine := $(shell uname -m)
LIBPATH = build/lib.linux-$(machine)-${python_version_major}.${python_version_minor}


help:
	@echo Available targets:
	@echo
	@echo '  pynauty       - build the pynauty extension module'
	@echo '  tests         - run all tests'
	@echo '  clean         - remove all python related temp files and dirs'
	@echo '  user-ins      - install pynauty into ~/.local/'
	@echo '  user-unins    - uninstall pynauty from ~/.local/'
	@echo '  virtenv-ins   - install pynauty into the active virtualenv'
	@echo '  virtenv-unins - uninstall pynauty from the active virtualenv'
	@echo '  dist          - create a source distribution'
	@echo '  docs          - build pyanauty documentation'
	@echo '  clean-docs    - remove pyanauty documentation'
	@echo '  nauty-objects - compile only nauty.o nautil.o naugraph.o schreier.o naurng.o'
	@echo '  nauty-progs   - build all nauty programs'
	@echo '  clean-nauty   - a "distclean" for nauty'
	@echo '  clobber       - clean + clean-nauty + clean-docs'
	@echo
	@echo Python version: ${python_version_full}
	@echo Machine type: ${machine}

pynauty: nauty-objects
	$(PYTHON) setup.py build

.PHONY: tests
tests: pynauty
	cd tests; PYTHONPATH="./../${LIBPATH}:$(PYTHONPATH)" $(PYTHON) test_autgrp.py
	cd tests; PYTHONPATH=".:../${LIBPATH}:$(PYTHONPATH)" $(PYTHON) test_isomorphic.py

virtenv-ins: pynauty
ifdef VIRTUAL_ENV
	$(PIP) install --upgrade .
else
	@echo ERROR: no VIRTUAL_ENV environment varaible found.
	@echo cannot install, aborting ...
	@exit 1
endif

virtenv-unins:
ifdef VIRTUAL_ENV
	$(PIP) uninstall pynauty
else
	@echo ERROR: no VIRTUAL_ENV environment varaible found.
	@echo cannot uninstall, aborting ...
	@exit 1
endif

user-ins: pynauty
	$(PIP) install --user --upgrade .

user-unins:
	$(PIP) uninstall pynauty

.PHONY: dist
dist: docs
	$(PYTHON) setup.py sdist

.PHONY: docs
docs:
	cd docs; make html

clean-docs:
	cd docs; make clean

clean:
	rm -fr build
	rm -fr dist
	rm -f MANIFEST
	rm -fr tests/{__pycache__,data_graphs.pyc}

# nauty stuff

GTOOLS = copyg listg labelg dretog amtog geng complg shortg showg NRswitchg \
  biplabg addedgeg deledgeg countg pickg genrang newedgeg catg genbg directg \
  multig planarg gentourng ranlabg runalltests subdivideg watercluster2 \
  linegraphg naucompare

nauty-config: $(NAUTY_DIR)/config.log
$(NAUTY_DIR)/config.log:
	cd $(NAUTY_DIR); ./configure CFLAGS='-O4 -fPIC'

nauty-objects: nauty-config
	cd $(NAUTY_DIR); make nauty.o nautil.o naugraph.o schreier.o naurng.o

nauty-programs: nauty-config
	cd $(NAUTY_DIR); make

clean-nauty:
	cd $(NAUTY_DIR); rm -f *.o dreadnaut ${GTOOLS} nauty.a nauty1.a
	cd $(NAUTY_DIR); rm -f makefile config.log config.status gtools.h naututil.h nauty.h

clobber: clean clean-nauty clean-docs

