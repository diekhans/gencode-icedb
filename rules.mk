# dependency file is generate as part of compile
${OBJDIR}/%.o: %.c
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -c -MM -MT $@ $< >$*.depend
	${CC} ${CFLAGS} -c -o $@ $<

(%.o): %.o
	@mkdir -p  $(dir $@)
	${ROOT}/make/addLib $@ $*.o

clean:
	rm -f ${PROGS} ${OBJS} ${LINKOBJS} ${DEPENDS}

savebak:
	savebak gencode-icedb *.h *.c Makefile

# don't fail on missing dependencies, they are first time the .o is generates
ifneq (${DEPENDS},)
-include ${DEPENDS}
endif

