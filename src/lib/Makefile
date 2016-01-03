ROOT = ../..
include ${ROOT}/config.mk

SRCS = intronMap.c \
	starSpliceJunction.c
OBJS = ${SRCS:%.c=${OBJDIR}/%.o}
DEPENDS = ${SRCS:%.c=%.depend}


all: ${OBJS} ${ICEDBLIB}(${OBJS})


include ${ROOT}/rules.mk

# Local Variables:
# mode: makefile-gmake
# End:
