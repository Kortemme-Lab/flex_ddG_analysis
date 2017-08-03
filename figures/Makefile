SHELL := /bin/bash

SOURCES = Makefile

all: fig-overview.pdf
fig-overview.pdf: fig-overview.dot $(SOURCES)
	dot -Tpdf -o fig-overview.pdf fig-overview.dot
clean:
	rm fig-overview.pdf
