###
# Generate variables from PROG_MODS.  Some of this is not pretty
###
BINPROGS = ${PROGS:%=${BINDIR}/%}

# define progName_PROG_OBJS
$(foreach prog,${PROGS},\
   $(eval ${prog}_PROG_OBJS=$(patsubst %,${OBJDIR}/%.o,${${prog}_PROG_MODS})))

# define progName_PROG_DEPENDS
$(foreach prog,${PROGS},\
   $(eval ${prog}_PROG_DEPENDS=$(patsubst %,.%.depend,${${prog}_PROG_MODS})))

# define LIB variables
LIB_DEPENDS = ${LIB_MODS:%=.%.depend}
LIB_OBJS = ${LIB_MODS:%=${OBJDIR}/%.o}

# define all objects and depends
OBJS = $(foreach prog,${PROGS},${${prog}_PROG_OBJS}) ${LIB_OBJS}
DEPENDS = $(foreach prog,${PROGS},${${prog}_PROG_DEPENDS}) ${LIB_DEPENDS}

# While compiling the .o's is straight forward, linking the programs is not
# because of not being able to generate the right dependencies.  So we link
# by recursion

all:: build

build:: ${PROGS:%=%_linkProg}

%_linkProg: ${OBJS}
	${MAKE} linkProg PROG=$*

ifneq (${PROG},)
# recursive call
linkProg: ${BINDIR}/${PROG}
${BINDIR}/${PROG}: $(${PROG}_PROG_OBJS) ${ICEDBLIB}
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -o $@ $^ ${LIBS}
endif

${OBJDIR}/%.o: %.c
	@mkdir -p $(dir $@)
	${CC} ${CFLAGS} -c -MM -MT $@ $< > .$*.depend
	${CC} ${CFLAGS} -c -o $@ $<

(%.o): %.o
	@mkdir -p  $(dir $@)
	${ROOT}/tools/addLib $@ $*.o

# don't fail on missing dependencies, they are generated the first time the .o
# is compiled
ifneq (${DEPENDS},)
-include ${DEPENDS}
endif

clean::
	rm -rf ${BINPROGS} ${OBJS} ${DEPENDS}
