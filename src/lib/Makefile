ROOT = ../..
include ${ROOT}/config.mk

SRCS =  intronMap.c \
	gencodeAttributes.c \
	starSpliceJunction.c
OBJS = ${SRCS:%.c=${OBJDIR}/%.o}
DEPENDS = ${SRCS:%.c=%.depend}


all: ${OBJS} ${ICEDBLIB}(${OBJS})

clean::
	rm -f ${ICEDBLIB}

include ${ROOT}/rules.mk

# Local Variables:
# mode: makefile-gmake
# End:
