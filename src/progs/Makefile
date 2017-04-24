ROOT = ../..

all: general.all rsl.all tsl.all

%.all:
	(cd $* && ${MAKE} all)

clean: general.all rsl.clean tsl.clean

%.clean:
	(cd $* && ${MAKE} clean)
