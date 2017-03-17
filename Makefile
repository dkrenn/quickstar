.PHONY: test

SAGE=sage

test: test-quickstar.py test-partitioner.py

test-%: %
	SAGE_PATH=$$PWD $(SAGE) -t --force-lib $<
