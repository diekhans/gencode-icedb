ROOT = ../../..

# force sort to be consistent
export LC_ALL=C


binDir = ${ROOT}/bin

test:: unitTests
	@echo "run mondoTests for heavy-duty tests"

unitTests:
	python classifyUnitTests.py


mondoTests: pslMondoTest gpMondotest

pslMondoTest: mkdirs
	./bin/mondoTester hg38 all_mrna psl

gpMondotest: mkdirs
	./bin/mondoTester hg38 wgEncodeGencodeCompV26 genepred


mkdirs:
	@mkdir -p output

clean::
	rm -rf output
