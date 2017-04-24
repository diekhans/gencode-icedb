ROOT = ../..

SUBDIRS = faUcscToGencode genePredIntrons

all: ${SUBDIRS:%=%.all}

%.all:
	(cd $* && ${MAKE} all)

clean: ${SUBDIRS:%=%.clean}

%.clean:
	(cd $* && ${MAKE} clean)
