ROOT = ../../..
include ${ROOT}/config.mk

LINKSRCS = intronMap.c \
	starSpliceJunction.c
PROGSRCS = spliceJunctionSupportSummarize.c
SRCS = ${LIBSRCS} ${PROGSRCS}
OBJS = ${SRCS:%.c=${OBJDIR}/%.o}
LINKOBJS = ${LINKSRCS:%.c=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.c=%.depend}
PROGS = ${BINDIR}/spliceJunctionSupportSummarize

all: ${OBJS} ${PROGS}

${BINDIR}/%: ${OBJDIR}/%.o ${LINKOBJS}
	@mkdir -p ${BINDIR}
	${CC} ${CFLAGS} -o $@ $< ${LINKOBJS} ${LIBS}

include ${ROOT}/rules.mk

# Local Variables:
# mode: makefile-gmake
# End:
