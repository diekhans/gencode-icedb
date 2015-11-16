
starGenerateGenome = ../bin/starGenerateGenome
starSpliceJunctionMap = ../bin/starSpliceJunctionMap

numThreads=2
test: smallTests

# tests with tiny region chr22:19175560-19178840 and reads from
# wgEncodeCaltechRnaSeqGm12878R1x75dAlignsRep1V2.bam
smallTests: smallGenomeTest smallMapSamTest

smallGenomeTest: mkdirs
	rm -rf output/$@.genome
	nice ${starGenerateGenome} --numThreads=${numThreads} --tmpDir=output input/small-chr22.fa input/small-gencode.gtf 75 output/$@.genome
	touch output/$@.done

smallMapSamTest: mkdirs output/smallGenomeTest.done
	nice ${starSpliceJunctionMap} --numThreads=${numThreads} --tmpDir=output output/smallGenomeTest.genome input/small-rnaseq.sam output/$@/$@
	diff expected/smallMapSamTest.outSJ.tab expected/smallMapSamTest.outSJ.tab


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
