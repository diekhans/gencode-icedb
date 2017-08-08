ROOT = ../../..

# force sort to be consistent
export LC_ALL=C

diff = diff -u

binDir = ${ROOT}/bin
gencodeDbLoad = ${binDir}/gencodeDbLoad

test:: unitTests gencodeDbLoadTest
	@echo "run mondoTests for heavy-duty tests"

unitTests:
	python classifyUnitTests.py

gencodeDbLoadTest: mkdirs
	rm -f output/$@.db
	${gencodeDbLoad} output/$@.db input/set1.gencodeAttrsV19.tsv input/set1.gencodeCompV19.gp
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_ann' >output/$@.gencode_ann.gp
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_attrs' >output/$@.gencode_attrs.tsv
	${diff} expected/$@.gencode_ann.gp output/$@.gencode_ann.gp
	${diff} expected/$@.gencode_attrs.tsv output/$@.gencode_attrs.tsv

mondoTests: pslMondoTest gpMondoTest gencodeDbMondoTest

pslMondoTest: mkdirs
	./bin/mondoTester hg38 all_mrna psl

gpMondoTest: mkdirs
	./bin/mondoTester hg38 wgEncodeGencodeCompV26 genepred

gencodeDbMondoTest: mkdirs
	rm -f output/$@.db
	${gencodeDbLoad} --hgdb=hg38 output/$@.db wgEncodeGencodeAttrsV26 wgEncodeGencodeCompV26 wgEncodeGencodePseudoGeneV26

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
