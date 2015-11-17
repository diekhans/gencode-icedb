
starGenerateGenome = ../bin/starGenerateGenome
starSpliceJunctionMap = ../bin/starSpliceJunctionMap
faUcscToGencode = ../bin/faUcscToGencode

.SECONDARY:

numThreads=2
test: smallTests faTests

# tests with tiny region chr22:19175560-19178840 and reads from
# wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep1V2.bam
smallTests: smallMapSamTest smallMapBamTest

# since generate genome take so long, run only once between make clean
output/smallGenomeTest.done: 
	rm -rf output/$@.genome
	@mkdir -p output
	nice ${starGenerateGenome} --numThreads=${numThreads} --tmpDir=output input/small-chr22.fa input/small-gencode.gtf 75 output/smallGenomeTest.genome
	touch output/smallGenomeTest.done

smallMapSamTest: mkdirs output/smallGenomeTest.done
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome input/small-rnaseq.sam output/smallMapSamTestSJ.out.tab
	diff expected/smallMapSamTestSJ.out.tab output/smallMapSamTestSJ.out.tab

smallMapBamTest: mkdirs output/smallGenomeTest.done
	samtools view -hb input/small-rnaseq.sam > output/$@.bam
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome output/$@.bam output/smallMapBamTestSJ.out.tab
	diff expected/smallMapSamTestSJ.out.tab output/smallMapBamTestSJ.out.tab

# tests of creating genome fasta with genomes names
faTests: ucscGencodeFaTest

ucscGencodeFaTest: mkdirs
	${faUcscToGencode} input/ucscGenome.fa input/ucsc2GencodeNames.tsv output/$@.fa
	diff expected/$@.fa output/$@.fa

mkdirs:
	@mkdir -p output

clean:
	rm -rf output
