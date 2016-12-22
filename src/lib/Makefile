ROOT = ../..

LIB_SRCS =  intronMap.c \
	gencodeAttributes.c \
	starSpliceJunction.c

include ${ROOT}/config.mk


all: ${LIB_OBJS} ${ICEDBLIB}(${LIB_OBJS})

clean::
	rm -f ${ICEDBLIB}

include ${ROOT}/rules.mk

# Local Variables:
# mode: makefile-gmake
# End:
