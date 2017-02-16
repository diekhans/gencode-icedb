ROOT = ../..
include ${ROOT}/config.mk

all::

test:: genbankEvidence.test classify.test

%.test:
	(cd $* && ${MAKE} test)

clean:: genbankEvidence.clean

%.clean:
	(cd $* && ${MAKE} clean)
