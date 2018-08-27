ROOT = ../..
include ${ROOT}/config.mk

faUcscToGencode = ${BINDIR}/faUcscToGencode
genePredIntrons = ${BINDIR}/genePredIntrons
ucscGencodeDbLoad = ${BINDIR}/ucscGencodeDbLoad
featureMondoTester = ./bin/featureMondoTester

.SECONDARY:

hsUcscGenomeDb = hg38
hsGencodeVer = V28
hsGencodeCompTbl = wgEncodeGencodeComp${hsGencodeVer}

test:   featureUnitTests \
	ucscGencodeDbLoadTest \
	faTests \
	genePredIntronsTests
	@echo "Note: run make mondoTest for heavy-duty tests"

##
# test of features
##
featureUnitTests:
	${PYTHON} featureUnitTests.py
 
##
# tests of loading gencode data into sqlite
##
ucscGencodeDbLoadTest: mkdirs
	rm -f output/$@.db
	${ucscGencodeDbLoad} --attrs=input/gencodeAttrsV28.tsv --transcriptSource=input/gencodeTranscriptSourceV28.tsv --transcriptionSupportLevel=input/gencodeTranscriptionSupportLevelV28.tsv  --tags=input/gencodeTagV28.tsv --genes=input/gencodeCompV28.gp output/$@.db
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
genePredIntronsTests: genePredIntronsSmallTest

genePredIntronsSmallTest: mkdirs
	faToTwoBit  input/small-chr22.fa output/$@.2bit
	${genePredIntrons} input/small-gencode.gp output/$@.2bit output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

mondoTest: pslMondoTest gpMondoTest gencodeDbMondoTest

pslMondoTest: mkdirs
	time ${featureMondoTester} ${hsUcscGenomeDb} all_mrna psl

gpMondoTest: mkdirs
	time ${featureMondoTester} ${hsUcscGenomeDb} ${hsGencodeCompTbl} genepred

gencodeDbMondoTest: mkdirs
	rm -f output/$@.db
	time ${ucscGencodeDbLoad} --hgdb=${hsUcscGenomeDb} --version=${hsGencodeVer} output/$@.db

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
