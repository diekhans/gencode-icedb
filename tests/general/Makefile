ROOT = ../..

binDir = ${ROOT}/bin

faUcscToGencode = ${binDir}/faUcscToGencode
genePredIntrons = ${binDir}/genePredIntrons
gencodeDbLoad = ${binDir}/gencodeDbLoad
featureMondoTester = ./bin/featureMondoTester

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
	${gencodeDbLoad} --attrs=input/set1.gencodeAttrsV19.tsv --transcriptSource=input/set1.gencodeTransSourceV19.tsv --transcriptionSupportLevel=input/set1.gencodeTransSupportLevelV19.tsv  --tags=input/set1.gencodeTagV19.tsv --genes=input/set1.gencodeCompV19.gp output/$@.db
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_ann' >output/$@.gencode_ann.gp
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_attrs' >output/$@.gencode_attrs.tsv
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_transcript_source' >output/$@.gencode_trans_src.tsv
	sqlite3 -header -batch output/$@.db 'SELECT * FROM gencode_transcription_support_level' >output/$@.gencode_trans_support_level.tsv
	${diff} expected/$@.gencode_ann.gp output/$@.gencode_ann.gp
	${diff} expected/$@.gencode_attrs.tsv output/$@.gencode_attrs.tsv
	${diff} expected/$@.gencode_trans_src.tsv output/$@.gencode_trans_src.tsv
	${diff} expected/$@.gencode_trans_support_level.tsv output/$@.gencode_trans_support_level.tsv


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
	time ${featureMondoTester} hg38 all_mrna psl

gpMondoTest: mkdirs
	time ${featureMondoTester} hg38 wgEncodeGencodeCompV26 genepred

gencodeDbMondoTest: mkdirs
	rm -f output/$@.db
	time ${gencodeDbLoad} --hgdb=hg38 --version=V26 output/$@.db

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
