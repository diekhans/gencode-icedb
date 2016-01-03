ROOT = ../../..
include ${ROOT}/config.mk

PROGSRCS = faUcscToGencode.c
SRCS = ${PROGSRCS}
OBJS = ${SRCS:%.c=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.c=%.depend}
PROGS = ${BINDIR}/faUcscToGencode

all: ${OBJS} ${PROGS}

${BINDIR}/%: ${OBJDIR}/%.o  ${ICEDBLIB}
	@mkdir -p ${BINDIR}
	${CC} ${CFLAGS} -o $@ $< ${LIBS}

include ${ROOT}/rules.mk

# Local Variables:
# mode: makefile-gmake
# End:
