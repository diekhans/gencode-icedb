ROOT = ..

SUBDIRS = faUcscToGencode genePredIntrons spliceJunctionCollectEvidence

all: ${SUBDIRS:%=%.all}

%.all:
	(cd $* && ${MAKE} all)

clean: ${SUBDIRS:%=%.clean}

%.clean:
	(cd $* && ${MAKE} clean)
