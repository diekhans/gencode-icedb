ROOT = ..
include ${ROOT}/config.mk

subdirs = general rsl tsl

all::

test:: ${subdirs:%=%.test}

%.test:
	(cd $* && ${MAKE} test)

mondoTest: general.mondoTest

%.mondoTest:
	(cd $* && ${MAKE} mondoTest)


clean:: ${subdirs:%=%.clean}

%.clean:
	(cd $* && ${MAKE} clean)

