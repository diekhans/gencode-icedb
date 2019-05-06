ROOT = ../..
include ${ROOT}/config.mk

.SECONDARY:

hsUcscGenomeDb = hg38
hsGencodeVer = V28
hsGencodeCompTbl = wgEncodeGencodeComp${hsGencodeVer}

assemblyRep = input/GCF_000001405.38_GRCh38.p12_assembly_report.txt

gencodeV28Db = output/gencodeV28.db
hg38TwoBit = /hive/data/genomes/hg38/hg38.2bit
hg38TestBams = output/ucsc-mrnaV28.bam output/hg38-mm10.transMap.bam output/ont-cls.bam


ensemblDb = homo_sapiens_core_92_38
grch38Fa = ../data/grch38/GRCh38.fa.gz

test:   featureUnitTests \
	ucscGencodeDbLoadTest \
	faTests \
	loadPslAlignsTest
	@echo "Note: run make mondoTest for heavy-duty tests"
	@echo "Note: run make ucscEnsemblCmp to compare import from UCSC database (via GTF) or direct from Ensembl database"

##
# test of features
##
unitTestCases = EvidenceTests GenePredAnnotationTests EnsemblDbAnnotationTests
featureUnitTests: ${unitTestCases:%=%_run}

%_run: ${hg38TestBams} ${gencodeV28Db}
	${PYTHON} featureUnitTests.py $*

 
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

loadPslAlignsTest: mkdirs
	${icedbLoadPslAligns} --table=$@ input/ucsc-mrnaV28.psl output/$@.db
	cut -f 10 input/ucsc-mrnaV28.psl | sort > output/$@.file.ids
	sqlite3 -batch -readonly output/$@.db 'SELECT qName FROM $@ ORDER BY qName' > output/$@.db.ids
	diff output/$@.file.ids output/$@.db.ids

mondoTest: pslMondoTest gpMondoTest gencodeDbMondoTest

pslMondoTest: mkdirs
	time ${featureMondoTester} ${hsUcscGenomeDb} all_mrna psl

gpMondoTest: mkdirs
	time ${featureMondoTester} ${hsUcscGenomeDb} ${hsGencodeCompTbl} genepred

gencodeDbMondoTest: mkdirs
	rm -f output/$@.db
	time ${ucscGencodeDbLoad} --hgdb=${hsUcscGenomeDb} --version=${hsGencodeVer} output/$@.db

ucscEnsemblCmp: ${gencodeV28Db} mkdirs
	${ucscEnsemblTransCmp} ${assemblyRep} ${gencodeV28Db} ${hg38TwoBit} ${ensemblDb} ${grch38Fa} output/$@.out

${gencodeV28Db}:
	@mkdir -p output
	rm -f $@ $@.tmp
	${ucscGencodeDbLoad} --hgdb=${hsUcscGenomeDb} --version=${hsGencodeVer} $@.tmp
	mv -f $@.tmp $@

output/%.bam: input/%.sam
	samtools sort -O BAM -o $@ $<
	samtools index -b $@

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
