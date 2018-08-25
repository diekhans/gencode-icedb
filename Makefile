ROOT = .
include ${ROOT}/config.mk

pyprogs = $(shell file -F $$'\t' bin/* tests/*/bin/* | awk '/Python script/{print $$1}')

all::
	(cd src && ${MAKE})

test:: all
	(cd tests && ${MAKE} test)

mondoTest::
	(cd tests && ${MAKE} mondoTest)

lint:
	flake8 tests lib/gencode_icedb ${pyprogs}

clean::
	(cd src && ${MAKE} clean)
	(cd tests && ${MAKE} clean)
	rm -rf  ${OBJDIR}
	find . -type f -name '*.pyc' -exec rm -f {} \;
	find . -depth -type d -name __pycache__ -exec rmdir {} \;
