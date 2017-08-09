ROOT = ../..

binDir = ${ROOT}/bin

rslEncodeDccQuery = ${binDir}/rslEncodeDccQuery
faUcscToGencode = ${binDir}/faUcscToGencode
genePredIntrons = ${binDir}/genePredIntrons

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

diff = diff -u

export LC_ALL=C

.SECONDARY:


test:   dataConvertTests \
	dbLoadTests \
	sjsupBuildTests \
	faTests \
	genePredIntronsTests \
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
sqldump = ${sqldumpcmd} 'select * from $(1)' > output/$@.$(1).tsv && ${diff} expected/$@.$(1).tsv output/$@.$(1).tsv

# test loading run-info from SRA
sraRunInfoDbLoad1Test: mkdirs
	rm -f output/$@.db
	${rslSraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	$(call sqldump,RunMetadata)

mappingMetadataDbLoadTest: mkdirs
	rm -f output/$@.db
	${rslSraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	${rslMappingMetadataDbLoad} output/$@.db input/set1.hs.mapping-params.tsv input/set1.hs.stardir.tsv
	$(call sqldump,RunMetadata)
	$(call sqldump,MappingParameters)
	$(call sqldump,MappingMetadata)

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
	for j in output/mkStarMergeTest.mergework/*.jobs ; do bash -e $$j ; done
	${diff} expected/$@.merged.sjsup output/$@.mergework/001.merged/000.sjsup
	bgzip -c output/mkStarMergeTest.mergework/001.merged/000.sjsup >output/$@.sjsup.gz
	tabix -0 -p bed output/$@.sjsup.gz

##
# Test of support
##
gencodeSupportTests: supportCollectTest

define gencodeSupportTestSetup
rm -rf output/$@
mkdir -p output/$@
${gencodeDbLoad} output/$@/gencode.db input/set1.hs.attrs.v26.tsv input/set1.hs.comp.v26.gp input/set1.hs.pseudo.v26.gp
${rslSraRunInfoDbLoad} output/$@/sup.db input/set1.hs.run-info.tsv
${rslMappingMetadataDbLoad} output/$@/sup.db input/set1.hs.mapping-params.tsv input/set1.hs.stardir.tsv
bgzip -c expected/mkStarMergeTest.merged.sjsup > output/$@/sup.sjsup.gz
tabix -0 -p bed output/$@/sup.sjsup.gz
endef

supportCollectTest: mkdirs
	${gencodeSupportTestSetup}
	${rslGencodeCollectSupport} output/$@/gencode.db output/$@/sup.db output/$@.supp.tsv ENSG00000142609.17 ENSG00000187730.8 ENSG00000226944.1


##
# tests of creating genome fasta with gencode names
##
faTests: ucscGencodeFaTest

ucscGencodeFaTest: mkdirs
	${faUcscToGencode} input/ucscGenome.fa input/ucsc2GencodeNames.tsv output/$@.fa
	diff expected/$@.fa output/$@.fa

##
# tests of program to get introns for genePreds
# FIXME: drop
##
genePredIntronsTests: genePredIntronsSet1Test

genePredIntronsSet1Test: mkdirs
	faToTwoBit  input/small-chr22.fa output/$@.2bit
	${genePredIntrons} input/small-gencode.gp output/$@.2bit output/$@.tsv
	diff expected/$@.tsv output/$@.tsv

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
