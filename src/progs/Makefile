ROOT = ../..
include ${ROOT}/config.mk

progs = ${BINDIR}/faUcscToGencode

all: ${progs}


${BINDIR}/%: %.c
	@mkdir -p ${BINDIR}
	${CC} ${CFLAGS} -o $@ $< ${LIBS}


clean:
	rm -f ${progs}
