ROOT = ../..

binDir = ${ROOT}/bin

rslEncodeDccQuery = ${binDir}/rslEncodeDccQuery

rslRnaSeqIntronEvidBed = ${binDir}/rslRnaSeqIntronEvidBed
rslRnaSeqIntronEvidStats = ${binDir}/rslRnaSeqIntronEvidStats

rslSraRunInfoFilter = ${binDir}/rslSraRunInfoFilter
rslSraRunInfoDbLoad = ${binDir}/rslSraRunInfoDbLoad
rslMappingMetadataDbLoad = ${binDir}/rslMappingMetadataDbLoad

rslStarSjOutSplit = ${binDir}/rslStarSjOutSplit
rslMkStarSjOutSplits = ${binDir}/rslMkStarSjOutSplits
rslMkStarSjSupMergeJobs = ${binDir}/rslMkStarSjSupMergeJobs

gencodeDbLoad = ${binDir}/gencodeDbLoad
rslGencodeCollectSupport = ${binDir}/rslGencodeCollectSupport
rslGencodeCollectSupportMkJobs = ${binDir}/rslGencodeCollectSupportMkJobs
rslGencodeCollectSupportFinishJobs = ${binDir}/rslGencodeCollectSupportFinishJobs

rslGencodeCollectNovel = ${binDir}/rslGencodeCollectNovel
rslGencodeCollectNovelMkJobs = ${binDir}/rslGencodeCollectNovelMkJobs
rslGencodeCollectNovelFinishJobs = ${binDir}/rslGencodeCollectNovelFinishJobs

genomeSeqsHg38 = /hive/data/genomes/hg38/hg38.2bit

diff = diff -u

# remove {make out exists ... } from job file so it can be executed as commands
jobsToCmds = sed -Ee 's/\{check out exists (.+)\}/\1/'

export LC_ALL=C

.SECONDARY:


test:   dataConvertTests \
	dbLoadTests \
	sjsupBuildTests \
	gencodeSupportTests \
	novelTests \
	encodeDccTests \
	intronEvidBedTests \
	intronEvidStatsTests

##
# conversion of data
##
dataConvertTests: sraRunInfoFilterStdTest \
		  sraRunInfoFilterFullTest \
		  sraRunInfoFilterSubsetTest \
		  sraRunInfoFilterMissingTest

sraRunInfoFilterStdTest: mkdirs
	${rslSraRunInfoFilter} input/runinfo.hs.csv output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

sraRunInfoFilterFullTest: mkdirs
	${rslSraRunInfoFilter} --full input/runinfo.hs.csv output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

sraRunInfoFilterSubsetTest: mkdirs
	${rslSraRunInfoFilter} --selectSubsetTsv=input/runinfo.hs.subset.tsv input/runinfo.hs.csv output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

# missing because of filter or just not there
sraRunInfoFilterMissingTest: mkdirs
	if ${rslSraRunInfoFilter} --subsetRequired --selectSubsetTsv=input/runinfo.hs.missing.tsv input/runinfo.hs.csv output/$@.tsv 2>output/$@.err ; then echo "expected failure" >&2; exit 1 ; fi
	${diff} expected/$@.err output/$@.err

##
# loading of database
##f
dbLoadTests: sraRunInfoDbLoad1Test \
		mappingMetadataDbLoadTest

sqldumpcmd = sqlite3 -header -batch output/$@.db

# call function with table as argument.
sqldumpdiff = ${sqldumpcmd} 'select * from $(1)' > output/$@.$(1).tsv && ${diff} expected/$@.$(1).tsv output/$@.$(1).tsv

# test loading run-info from SRA
sraRunInfoDbLoad1Test: mkdirs
	rm -f output/$@.db
	${rslSraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	$(call sqldumpdiff,RunMetadata)

mappingMetadataDbLoadTest: mkdirs
	rm -f output/$@.db
	${rslSraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	${rslMappingMetadataDbLoad} output/$@.db input/set1.hs.mapping-params.tsv input/set1.hs.stardir.tsv
	$(call sqldumpdiff,RunMetadata)
	$(call sqldumpdiff,MappingParameters)
	$(call sqldumpdiff,MappingMetadata)

##
# build of tabix database
##
sjsupBuildTests: sjoutSplitTest \
		sjoutSplitSomeTest \
		mkStarSplitsTest \
		mkStarMergeTest

# md5 is saved to avoid having lots of expected files
sjoutSplitTest: mkdirs
	rm -rf output/$@.split
	${rslStarSjOutSplit} input/set1.hs.stardir.tsv output/$@.split
	openssl md5 output/$@.split/*/*.sjsup > output/$@.md5
	${diff} expected/$@.md5 output/$@.md5

sjoutSplitSomeTest: mkdirs
	rm -rf output/$@.split
	${rslStarSjOutSplit} -startDirRec=0 -endDirRec=1 input/set1.hs.stardir.tsv output/$@.split
	openssl md5 output/$@.split/*/*.sjsup > output/$@.md5
	${diff} expected/$@.md5 output/$@.md5

mkStarSplitsTest: mkdirs
	rm -rf output/$@.mergework
	${rslMkStarSjOutSplits} --numprocs=2 input/set1.hs.stardir.tsv output/$@.mergework
	openssl md5 output/$@.mergework/byChrom/*/*.sjsup > output/$@.md5
	${diff} expected/$@.md5 output/$@.md5

# path in job file varies, just check for count (and make wc output portable)
mkStarMergeTest: mkdirs
	rm -rf output/$@.mergework
	${rslMkStarSjOutSplits} --numprocs=2 input/set1.hs.stardir.tsv output/$@.mergework
	${rslMkStarSjSupMergeJobs} --filesPerJob=5 output/$@.mergework
	wc -l output/$@.mergework/*.jobs | tawko '{print $$2,$$1}' > output/$@.jobs.cnt
	${diff} expected/$@.jobs.cnt output/$@.jobs.cnt
	${diff} expected/$@.run-batches output/$@.mergework/run-batches
	for j in output/mkStarMergeTest.mergework/*.jobs ; do ${jobsToCmds} $$j | bash -e ; done
	${diff} expected/$@.merged.sjsup output/$@.mergework/001.merged/000.sjsup
	bgzip -c output/mkStarMergeTest.mergework/001.merged/000.sjsup >output/$@.sjsup.gz
	tabix -0 -p bed output/$@.sjsup.gz

##
# Test of support
##
gencodeSupportTests: supportCollectAllTest \
			supportCollectGenesTest \
			supportCollectTransTest \
			supportCollectMkJobsTest

define gencodeSupportTestSetup
rm -rf output/$@
mkdir -p output/$@
${gencodeDbLoad} output/$@/gencode.db input/set1.hs.attrs.v26.tsv input/set1.hs.trans-source.v26.tab input/set1.hs.comp.v26.gp input/set1.hs.pseudo.v26.gp
${rslSraRunInfoDbLoad} output/$@/rnaseq.db input/set1.hs.run-info.tsv
${rslMappingMetadataDbLoad} output/$@/rnaseq.db input/set1.hs.mapping-params.tsv input/set1.hs.stardir.tsv
bgzip -c expected/mkStarMergeTest.merged.sjsup > output/$@/rnaseq.sjsup.gz
tabix -0 -p bed output/$@/rnaseq.sjsup.gz
endef

# uncomment to test mock sequences
#forceMockGenomeSeqs = --forceMockGenomeSeqs

# uncomment to update mock in supportCollectAllTest
#updateMockGenomeSeqs = --updateMockGenomeSeqs

genomeSeqsSpecs = --genomeSeqs=${genomeSeqsHg38} --mockGenomeSeqs=input/set1.hs.mock-genome.tsv


supportCollectAllTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectSupport} ${genomeSeqsSpecs} ${updateMockGenomeSeqs} ${forceMockGenomeSeqs} output/$@/gencode.db output/$@/rnaseq.db output/$@.supp.tsv
	${diff} expected/$@.supp.tsv output/$@.supp.tsv

supportCollectGenesTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectSupport} ${genomeSeqsSpecs} ${forceMockGenomeSeqs} output/$@/gencode.db output/$@/rnaseq.db output/$@.supp.tsv ENSG00000235169.7
	${diff} expected/$@.supp.tsv output/$@.supp.tsv

supportCollectTransTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectSupport} ${genomeSeqsSpecs} ${forceMockGenomeSeqs} output/$@/gencode.db output/$@/rnaseq.db output/$@.supp.tsv ENST00000444870.6 ENST00000462356.5
	${diff} expected/$@.supp.tsv output/$@.supp.tsv

supportCollectMkJobsTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectSupportMkJobs} --maxgenes=10 ${genomeSeqsSpecs} output/$@/gencode.db output/$@/rnaseq.db output/$@.work
	sed -Ee 's#^.*(/bin/rslGencodeCollectSupport)#??\1#' output/$@.work/batch.jobs > output/$@.work/batch.jobs.cleaned
	${diff} expected/$@.batch.jobs.cleaned output/$@.work/batch.jobs.cleaned
	${jobsToCmds} output/$@.work/batch.jobs | bash -e  
	ls -1 output/$@.work/results/ > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${rslGencodeCollectSupportFinishJobs} output/$@.work output/$@.db
	$(call sqldumpdiff,GencodeSupport)

##
# test of novel
##
novelTests: novelRange1Test \
	novelMkJobsTest

novelRange1Test: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectNovel} ${genomeSeqsSpecs} output/$@/gencode.db output/$@/rnaseq.db output/$@.novel.tsv chr1:2193910-7731544
	${diff} expected/$@.novel.tsv output/$@.novel.tsv

# make sure putative introns start in window
novelSplitTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectNovel} ${genomeSeqsSpecs} output/$@/gencode.db output/$@/rnaseq.db output/$@.win1.novel.tsv chr22:0-10940708
	${rslGencodeCollectNovel} ${genomeSeqsSpecs} output/$@/gencode.db output/$@/rnaseq.db output/$@.win2.novel.tsv chr22:10940708-20000000
	${diff} expected/$@.win1.novel.tsv output/$@.win1.novel.tsv
	${diff} expected/$@.win2.novel.tsv output/$@.win2.novel.tsv

novelMkJobsTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectNovelMkJobs} --testRange=chr22:10640100-17601091 --windowSize=1000000 ${genomeSeqsSpecs} output/$@/gencode.db output/$@/rnaseq.db output/$@.work
	sed -Ee 's#^.*(/bin/rslGencodeCollectNovel)#??\1#' output/$@.work/batch.jobs > output/$@.work/batch.jobs.cleaned
	${diff} expected/$@.batch.jobs.cleaned output/$@.work/batch.jobs.cleaned
	${jobsToCmds} output/$@.work/batch.jobs | bash -e  
	ls -1 output/$@.work/results/ > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${rslGencodeCollectNovelFinishJobs} output/$@.work output/$@.db
	$(call sqldumpdiff,GencodeNovel)


##
# data from encode3 DCC
##
encodeDccTests: encodeDcc10Tests

# just test aginst live site, it was too much of a pain to edit
# down the huge response to use it as fake input.  Can't verify
# output.
encodeDcc10Tests: mkdirs
	${rslEncodeDccQuery} "Mus musculus" --limit=10 output/encodeDcc10Tests.tsv output/encodeDcc10Tests.urls

##
# rslRnaSeqIntronEvidBed
##
intronEvidBedTests: intronEvidBed1Test

intronEvidBed1Test: mkdirs
	${rslRnaSeqIntronEvidBed} input/hs.chr22-2k.evid.tsv output/$@.bed

##
# rslRnaSeqIntronEvidStats
##
intronEvidStatsTests: intronEvidStatsAllTest intronEvidStatsMin30 \
	intronEvidStatsConsensus

intronEvidStatsAllTest: mkdirs
	${rslRnaSeqIntronEvidStats} input/hs.chr22-2k.evid.tsv output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv

intronEvidStatsMin30: mkdirs
	${rslRnaSeqIntronEvidStats} --minIntronSize=30 input/hs.chr22-2k.evid.tsv output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv

intronEvidStatsConsensus: mkdirs
	${rslRnaSeqIntronEvidStats} --spliceJuncCat=consensus input/hs.chr22-2k.evid.tsv output/$@.stats.tsv
	diff expected/$@.stats.tsv output/$@.stats.tsv


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
