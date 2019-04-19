ROOT = ../..

SUBDIRS = faUcscToGencode

# not used: genePredIntrons

all: ${SUBDIRS:%=%.all}

%.all:
	(cd $* && ${MAKE} all)

clean: ${SUBDIRS:%=%.clean}

%.clean:
	(cd $* && ${MAKE} clean)
