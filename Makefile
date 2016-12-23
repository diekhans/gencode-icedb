ROOT = .
include ${ROOT}/config.mk

all::
	(cd src && ${MAKE})

test::
	(cd tests && ${MAKE} test)

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf ${BINDIR}/*.dSYM ${OBJS} lib/gencode_icedb/*.pyc objs
