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
dataConvertTests: sraRunInfoFilterTest \
		  sraRunInfoFilterSubsetTest \
		  sraRunInfoFilterMissingTest

sraRunInfoFilterTest: mkdirs
	${sraRunInfoFilter} input/runinfo.hs.csv output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

sraRunInfoFilterSubsetTest: mkdirs
	${sraRunInfoFilter} --selectSubsetTsv=input/runinfo.hs.subset.tsv input/runinfo.hs.csv output/$@.tsv
	${diff} expected/$@.tsv output/$@.tsv

# missing because of filter or just not there
sraRunInfoFilterMissingTest: mkdirs
	if ${sraRunInfoFilter} --subsetRequired --selectSubsetTsv=input/runinfo.hs.missing.tsv input/runinfo.hs.csv output/$@.tsv 2>output/$@.err ; then echo "expected failure" >&2; exit 1 ; fi
	${diff} expected/$@.err output/$@.err

##
# loading of database
##f
dbLoadTests: sraRunInfoDbLoad1Test \
		mappingMetadataDbLoadTest

# test loading run-info from SRA
sraRunInfoDbLoad1Test: mkdirs
	rm -f output/$@.db
	${sraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	$(call sqldumpdiff,run_metadata)

mappingMetadataDbLoadTest: mkdirs
	rm -f output/$@.db
	${sraRunInfoDbLoad} output/$@.db input/set1.hs.run-info.tsv
	${rslMappingMetadataDbLoad} output/$@.db input/set1.hs.mapping-params.tsv input/set1.hs.stardir.tsv
	$(call sqldumpdiff,run_metadata)
	$(call sqldumpdiff,mapping_parameters)
	$(call sqldumpdiff,mapping_metadata)

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
