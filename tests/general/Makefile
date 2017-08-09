ROOT = ../..

binDir = ${ROOT}/bin

faUcscToGencode = ${binDir}/faUcscToGencode
genePredIntrons = ${binDir}/genePredIntrons
gencodeDbLoad = ${binDir}/gencodeDbLoad

diff = diff -u

export LC_ALL=C

.SECONDARY:


test:   featureUnitTests \
	gencodeDbLoadTest \
	faTests \
	genePredIntronsTests
	@echo "Note: run make mondoTest for heavy-duty tests"

##
# test of features
##
featureUnitTests:
	python featureUnitTests.py

##
# tests of loading gencode data into sqlite
##
gencodeDbLoadTest: mkdirs
	rm -f output/$@.db
	${gencodeDbLoad} output/$@.db input/set1.gencodeAttrsV19.tsv input/set1.gencodeCompV19.gp
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_ann' >output/$@.gencode_ann.gp
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_attrs' >output/$@.gencode_attrs.tsv
	${diff} expected/$@.gencode_ann.gp output/$@.gencode_ann.gp
	${diff} expected/$@.gencode_attrs.tsv output/$@.gencode_attrs.tsv


##
# tests of creating genome fasta with gencode names
##
faTests: ucscGencodeFaTest

ucscGencodeFaTest: mkdirs
	${faUcscToGencode} input/ucscGenome.fa input/ucsc2GencodeNames.tsv output/$@.fa
	diff expected/$@.fa output/$@.fa

##
# tests of program to get introns for genePreds
##
genePredIntronsTests: genePredIntronsSet1Test

genePredIntronsSet1Test: mkdirs
	faToTwoBit  input/small-chr22.fa output/$@.2bit
	${genePredIntrons} input/small-gencode.gp output/$@.2bit output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

mondoTest: pslMondoTest gpMondoTest gencodeDbMondoTest

pslMondoTest: mkdirs
	./bin/featureMondoTester hg38 all_mrna psl

gpMondoTest: mkdirs
	./bin/featureMondoTester hg38 wgEncodeGencodeCompV26 genepred

gencodeDbMondoTest: mkdirs
	rm -f output/$@.db
	${gencodeDbLoad} --hgdb=hg38 output/$@.db wgEncodeGencodeAttrsV26 wgEncodeGencodeCompV26 wgEncodeGencodePseudoGeneV26

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
