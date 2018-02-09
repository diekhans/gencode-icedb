ROOT = ../..
include ${ROOT}/config.mk

subdirs = evidence classify

all::

test:: ${subdirs:%=%.test}

%.test:
	(cd $* && ${MAKE} test)

clean:: ${subdirs:%=%.clean}

%.clean:
	(cd $* && ${MAKE} clean)
