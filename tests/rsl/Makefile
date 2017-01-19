ROOT = ../..

binDir = ${ROOT}/bin
starGenerateGenome = ${binDir}/starGenerateGenome
starSpliceJunctionMap = ${binDir}/starSpliceJunctionMap
faUcscToGencode = ${binDir}/faUcscToGencode
estimateReadLength = ${binDir}/estimateReadLength
rnaSeqDataRegister = ${binDir}/rnaSeqDataRegister
spliceJunctionCollectEvidence = ${binDir}/spliceJunctionCollectEvidence
genePredIntrons = ${binDir}/genePredIntrons
encodeDccQuery = ${binDir}/encodeDccQuery
rnaSeqIntronSupport = ${binDir}/rnaSeqIntronSupport
rnaSeqIntronEvidBed = ${binDir}/rnaSeqIntronEvidBed
rnaSeqIntronEvidStats = ${binDir}/rnaSeqIntronEvidStats
sjCollectEvidence = ${binDir}/sjCollectEvidence

.SECONDARY:

numThreads=2
test: smallTests readLenTests faTests rnaSeqDataRegisterTest genePredIntronsTest \
	spliceJuncEvidenceTests sjCollectEvidenceTests encodeDccTests intronSupportTests \
	intronEvidBedTests intronEvidStatsTests

##
# small test data:
##
smallTests: smallGenomeTest smallGenomeCompressedTest smallMapSamTest smallMapBamTest smallMapPairFastqTest

smallGenomeTest: mkdirs
	rm -rf output/$@.genome
	nice ${starGenerateGenome} --numThreads=${numThreads} --tmpDir=output input/small-chr22.fa input/small-gencode.gtf 75 output/$@.genome

smallGenomeCompressedTest: mkdirs
	rm -rf output/$@.genome
	gzip -c input/small-gencode.gtf > output/$@.gtf.gz
	gzip -c input/small-chr22.fa > output/$@.fa.gz
	nice ${starGenerateGenome} --numThreads=${numThreads} --tmpDir=output output/$@.fa.gz output/$@.gtf.gz 75 output/$@.genome


# split test reads file into a pair of paried read fastq files
output/smallFastqPair.done:
	@mkdir -p output
	samtools sort -n --output-fmt sam input/small-rnaseq.sam  |samtools fastq -1 output/smallFastqPair_1.fastq -2 output/smallFastqPair_2.fastq -s output/smallFastqPair_0.fastq /dev/stdin
	touch $@

smallMapSamTest: mkdirs smallGenomeTest
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome input/small-rnaseq.sam output/$@.sj.tab
	diff expected/$@.sj.tab output/$@.sj.tab

smallMapBamTest: mkdirs smallGenomeTest
	samtools view -hb input/small-rnaseq.sam > output/$@.bam
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome output/$@.bam output/$@.sj.tab
	diff expected/smallMapSamTest.sj.tab output/$@.sj.tab

smallMapPairFastqTest: mkdirs smallGenomeTest output/smallFastqPair.done
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome output/smallFastqPair_1.fastq --readsFile2=output/smallFastqPair_2.fastq output/$@.sj.tab
	diff expected/$@.sj.tab output/$@.sj.tab
##
# rest obtaining average read lengths
##
readLenTests: readLenAllTest readLenSomeTest

readLenAllTest: mkdirs
	${estimateReadLength} input/small-rnaseq.sam output/$@.avelen
	diff expected/$@.avelen output/$@.avelen

readLenSomeTest: mkdirs
	${estimateReadLength} --sampleSize=5 input/small-rnaseq.sam > output/$@.avelen
	diff expected/$@.avelen output/$@.avelen
##
# tests of creating genome fasta with genomes names
##
faTests: ucscGencodeFaTest

ucscGencodeFaTest: mkdirs
	${faUcscToGencode} input/ucscGenome.fa input/ucsc2GencodeNames.tsv output/$@.fa
	diff expected/$@.fa output/$@.fa
##
# setup data directory and files for tests using data.db.
##
define dataTreeSetup
	@rm -rf output/$@
	@mkdir -p output/$@/data/human/rnaSeq/smallSet output/$@/data/human/rnaSeq/otherSet
	cp input/small-rnaseq.sam output/$@/data/human/rnaSeq/smallSet/small.sam
	samtools view -hb input/small-rnaseq.sam > output/$@/data/human/rnaSeq/smallSet/small.bam
	ln -f output/smallFastqPair_1.fastq output/smallFastqPair_2.fastq output/$@/data/human/rnaSeq/smallSet/
	ln -f output/smallFastqPair_1.fastq output/smallFastqPair_2.fastq output/$@/data/human/rnaSeq/otherSet/
	${rnaSeqDataRegister} --rootDir=output/$@ \
		--readsFileUrl="http://fred.org/data/input/small-rnaseq.sam" \
		--description="small RNA-Seq" \
		--tissue="micro brain" \
		human smallSet small1A small.sam
	${rnaSeqDataRegister} --rootDir=output/$@ \
		--readsFileUrl="http://fred.org/data/input/smallish-rnaseq.bam" \
		--description="not so small RNA-Seq" \
		--tissue="macro brain" \
		human smallSet small1B small.bam
	${rnaSeqDataRegister} --rootDir=output/$@ --skipExisting \
		human smallSet small1B small.bam
	${rnaSeqDataRegister} --rootDir=output/$@ --skipExisting \
		--readsFileUrl="http://fred.org/data/input/smallFastqPair_1.fastq" \
		--readsFile2Url="http://fred.org/data/input/smallFastqPair_2.fastq" \
		--description="paired RNA-Seq" \
		--tissue="micro macro brain" \
		human smallSet small1-2 smallFastqPair_1.fastq \
		--readsFile2=smallFastqPair_2.fastq
endef

# need to construct data test in the expected locations
rnaSeqDataRegisterTest: mkdirs output/smallFastqPair.done
	${dataTreeSetup}
	echo 'select * from rnaseqdata;' | sqlite3 -bail output/$@/data/data.db > output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

##
# tests of program to get introns for genePreds
##
genePredIntronsTest: genePredIntronsSet1Test

genePredIntronsSet1Test: mkdirs
	faToTwoBit  input/small-chr22.fa output/$@.2bit
	${genePredIntrons} input/small-gencode.gp output/$@.2bit output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

##
# tests splice junction evidence collection
##
spliceJuncEvidenceTests: spliceJuncEvidBasicTest spliceJuncEvidMinOverhangTest \
		spliceJuncEvidNoSupportTest

spliceJuncEvidBasicTest: mkdirs
	${spliceJunctionCollectEvidence} input/small-gencode.gp expected/genePredIntronsSet1Test.tsv input/small-sj.lst output/$@.evid.tsv
	diff expected/$@.evid.tsv output/$@.evid.tsv

spliceJuncEvidMinOverhangTest: mkdirs
	${spliceJunctionCollectEvidence} -minOverhang=20 input/small-gencode.gp expected/genePredIntronsSet1Test.tsv input/small-sj.lst output/$@.evid.tsv
	diff expected/$@.evid.tsv output/$@.evid.tsv

# edit to remove some evidence to check report on unsupport introns
spliceJuncEvidNoSupportTest: mkdirs output/smallMapPairFastqTest-edit.sj.tab
	${spliceJunctionCollectEvidence} input/small-gencode.gp expected/genePredIntronsSet1Test.tsv input/small-sj-edit.lst output/$@.evid.tsv
	diff expected/$@.evid.tsv output/$@.evid.tsv

##
# tests of program that runs splice junction evidence collection based on
# database
##
sjCollectEvidenceTests: sjCollectEvidMinOverhangTest

sjCollectEvidMinOverhangSjDir = output/sjCollectEvidMinOverhangTest/analysis/grch38/smallSet/rnaSeqEvid/smallSet
sjCollectEvidMinOverhangTest: mkdirs output/smallFastqPair.done
	${dataTreeSetup}
	@mkdir -p ${sjCollectEvidMinOverhangSjDir}
	cp expected/smallMapPairFastqTest.sj.tab ${sjCollectEvidMinOverhangSjDir}/small1A.sj.tab
	cp expected/smallMapSamTest.sj.tab ${sjCollectEvidMinOverhangSjDir}/small1B.sj.tab
	cp expected/smallMapPairFastqTest.sj.tab ${sjCollectEvidMinOverhangSjDir}/small1-2.sj.tab
	${sjCollectEvidence} --rootDir=output/$@ --minOverhang=20 human grch38 smallSet input/small-gencode.gp expected/genePredIntronsSet1Test.tsv
	diff expected/$@.evid.tsv output/sjCollectEvidMinOverhangTest/analysis/grch38/smallSet/rnaSeqSupport/intron.evid.tsv


# filtered files that has some introns unsuppported
output/smallMapPairFastqTest-edit.sj.tab: expected/smallMapPairFastqTest.sj.tab
	@mkdir -p output
	awk 'BEGIN{IFS=FS="\\t"} $$2<2000' expected/smallMapPairFastqTest.sj.tab > output/smallMapPairFastqTest-edit.sj.tab

##
encodeDccTests: encodeDcc10Tests

# just test aginst live site, it was too much of a pain to edit
# down the huge response to use it as fake input.  Can't verify
# output.
encodeDcc10Tests: mkdirs
	${encodeDccQuery} "Mus musculus" --limit=10 output/encodeDcc10Tests.tsv output/encodeDcc10Tests.urls

##
# rnaSeqIntronSupport tests
##
intronSupportTests: intronSupport1Test

intronSupport1Test: mkdirs
	${rnaSeqIntronSupport} input/small-gencode.gp input/small-gencode.attrs \
		expected/spliceJuncEvidNoSupportTest.evid.tsv

##
# rnaSeqIntronEvidBed
##
intronEvidBedTests: intronEvidBed1Test

intronEvidBed1Test: mkdirs
	${rnaSeqIntronEvidBed} input/hs.chr22-2k.evid.tsv output/$@.bed

##
# rnaSeqIntronEvidStats
##
intronEvidStatsTests: intronEvidStatsAllTest intronEvidStatsMin30 \
	intronEvidStatsConsensus

intronEvidStatsAllTest: mkdirs
	${rnaSeqIntronEvidStats} input/hs.chr22-2k.evid.tsv output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv

intronEvidStatsMin30: mkdirs
	${rnaSeqIntronEvidStats} --minIntronSize=30 input/hs.chr22-2k.evid.tsv output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv

intronEvidStatsConsensus: mkdirs
	${rnaSeqIntronEvidStats} --spliceJuncCat=consensus input/hs.chr22-2k.evid.tsv output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
