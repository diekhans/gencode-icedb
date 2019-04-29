ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/rsl/rslTestDefs.mk



test:   dataConvertTests \
	dbLoadTests \
	sjsupBuildTests \
	encodeDccTests

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

# test loading run-info from SRA
sraRunInfoDbLoad1Test: mkdirs
	rm -f output/$@.db
	${rslSraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	$(call sqldumpdiff,run_metadata)

mappingMetadataDbLoadTest: mkdirs
	rm -f output/$@.db
	${rslSraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	${rslMappingMetadataDbLoad} output/$@.db input/set1.hs.mapping-params.tsv input/set1.hs.stardir.tsv
	$(call sqldumpdiff,run_metadata)
	$(call sqldumpdiff,mapping_parameters)
	$(call sqldumpdiff,mapping_metadata)

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
# batches must be run in order, hence muliple jobsToMake
mkStarMergeTest: mkdirs
	rm -rf output/$@.mergework
	${rslMkStarSjOutSplits} --numprocs=2 input/set1.hs.stardir.tsv output/$@.mergework
	${rslMkStarSjSupMergeJobs} --filesPerJob=5 output/$@.mergework
	wc -l output/$@.mergework/*.jobs | tawko '{print $$2,$$1}' > output/$@.jobs.cnt
	${diff} expected/$@.jobs.cnt output/$@.jobs.cnt
	${diff} expected/$@.run-batches output/$@.mergework/run-batches
	for j in output/mkStarMergeTest.mergework/*.jobs ; do ${jobsToMake} $$j output/$@.mk && ${MAKE} -f output/$@.mk ; done
	${diff} expected/$@.merged.sjsup output/$@.mergework/001.merged/000.sjsup
	bgzip -c output/mkStarMergeTest.mergework/001.merged/000.sjsup >output/$@.sjsup.gz
	tabix -0 -p bed output/$@.sjsup.gz

##
# data from encode3 DCC
##
encodeDccTests: encodeDcc10Tests

# just test aginst live site, it was too much of a pain to edit
# down the huge response to use it as fake input.  Can't verify
# output.
encodeDcc10Tests: mkdirs
	${rslEncodeDccQuery} "Mus musculus" --limit=10 output/encodeDcc10Tests.tsv output/encodeDcc10Tests.urls


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
