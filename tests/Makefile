
starGenerateGenome = ../bin/starGenerateGenome
starSpliceJunctionMap = ../bin/starSpliceJunctionMap
faUcscToGencode = ../bin/faUcscToGencode
estimateReadLength = ../bin/estimateReadLength
rnaSeqDataRegister = ../bin/rnaSeqDataRegister

.SECONDARY:

numThreads=2
test: smallTests readLenTests faTests


# small test data:

smallTests: smallMapSamTest smallMapBamTest smallMapPairFastqTest

# since generate genome take so long, run only once between make clean
output/smallGenomeTest.done: 
	rm -rf output/smallGenomeTest.genome
	@mkdir -p output
	nice ${starGenerateGenome} --numThreads=${numThreads} --tmpDir=output input/small-chr22.fa input/small-gencode.gtf 75 output/smallGenomeTest.genome
	touch $@

# split test reads file into a pair of paried read fastq files
output/smallFastqPair.done:
	@mkdir -p output
	samtools sort -n --output-fmt sam input/small-rnaseq.sam  |samtools fastq -1 output/smallFastqPair_1.fastq -2 output/smallFastqPair_2.fastq -s output/smallFastqPair_0.fastq /dev/stdin
	touch $@

smallMapSamTest: mkdirs output/smallGenomeTest.done
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome input/small-rnaseq.sam output/$@SJ.out.tab
	diff expected/$@SJ.out.tab output/$@SJ.out.tab

smallMapBamTest: mkdirs output/smallGenomeTest.done
	samtools view -hb input/small-rnaseq.sam > output/$@.bam
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome output/$@.bam output/$@SJ.out.tab
	diff expected/smallMapSamTestSJ.out.tab output/$@SJ.out.tab

smallMapPairFastqTest: mkdirs output/smallGenomeTest.done output/smallFastqPair.done
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome output/smallFastqPair_1.fastq --readsFile2=output/smallFastqPair_2.fastq output/$@SJ.out.tab
	diff expected/$@SJ.out.tab output/$@SJ.out.tab

# rest obtaining average read lengths
readLenTests: readLenAllTest readLenSomeTest

readLenAllTest: mkdirs
	${estimateReadLength} input/small-rnaseq.sam output/$@.avelen
	diff expected/$@.avelen output/$@.avelen

readLenSomeTest: mkdirs
	${estimateReadLength} --sampleSize=5 input/small-rnaseq.sam > output/$@.avelen
	diff expected/$@.avelen output/$@.avelen

# tests of creating genome fasta with genomes names
faTests: ucscGencodeFaTest

ucscGencodeFaTest: mkdirs
	${faUcscToGencode} input/ucscGenome.fa input/ucsc2GencodeNames.tsv output/$@.fa
	diff expected/$@.fa output/$@.fa

rnaSeqDataRegisterTest: mkdirs
	@rm -f output/$@.db
	samtools view -hb input/small-rnaseq.sam > output/$@.bam
	${rnaSeqDataRegister} --dataRootDir=. \
		--readsFileUrl="http://fred.org/data/input/small-rnaseq.sam" \
		--description="small RNA-Seq" \
		--tissue="micro brain" \
		small1A "human" input/small-rnaseq.sam output/$@.db
	${rnaSeqDataRegister} --dataRootDir=. \
		--readsFileUrl="http://fred.org/data/input/smallish-rnaseq.bam" \
		--description="not so small RNA-Seq" \
		--tissue="macro brain" \
		small1B "human" output/$@.bam output/$@.db
	echo 'select * from rnaseqdata;' | sqlite3 -bail output/$@.db > output/$@.tsv
	cut -f 1-2,4- output/$@.tsv >output/$@.notime.tsv
	diff expected/$@.notime.tsv output/$@.notime.tsv

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
