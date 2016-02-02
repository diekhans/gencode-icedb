ROOT = ../../..
include ${ROOT}/config.mk

LINKSRCS = intronCounts.c
PROGSRCS = spliceJunctionCollectEvidence.c
SRCS = ${LINKSRCS} ${PROGSRCS}
OBJS = ${SRCS:%.c=${OBJDIR}/%.o}
LINKOBJS = ${LINKSRCS:%.c=${OBJDIR}/%.o}
DEPENDS =  ${SRCS:%.c=%.depend}
PROGS = ${BINDIR}/spliceJunctionCollectEvidence

all: ${OBJS} ${PROGS}

${BINDIR}/%: ${OBJDIR}/%.o ${LINKOBJS} ${ICEDBLIB}
	@mkdir -p ${BINDIR}
	${CC} ${CFLAGS} -o $@ $< ${LINKOBJS} ${LIBS}

include ${ROOT}/rules.mk

# Local Variables:
# mode: makefile-gmake
# End:
