SHELL := /bin/bash

SOURCES = Makefile

all: fig1.pdf
fig1.pdf: fig1.dot $(SOURCES)
	dot -Tpdf -o fig1.pdf fig1.dot
clean:
	rm fig1.pdf
