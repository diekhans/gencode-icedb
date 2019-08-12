ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/rsl/rslTestDefs.mk

test: arrayExpressTests


test2:   dataConvertTests \
	dbLoadTests \
	sjsupBuildTests \
	encodeDccTests

##
# data from array express rnaseqer
##
arrayExpressTests: arrayExpressSourceTest

arrayExpressSourceTest: mkdirs
	rm -f output/$@.db
	${rslArrayExpressSourceLoad} \
		--analysis=input/ERR2019438/ERR2019438.rnaseqer.json \
		--analysis=input/ERR3301312/ERR3301312.rnaseqer.json \
		--analysis=input/SRR3555865/SRR3555865.rnaseqer.json \
		output/$@.db
	$(call sqldumpcmd)  'select * from rsl_evidence_source' | tawk 'NR>1{$$2 = $$3 = $$4 = $$9 = "XXX"}{print $$0}' > output/$@.rsl_evidence_source.tsv
	${diff} expected/$@.rsl_evidence_source.tsv output/$@.rsl_evidence_source.tsv

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
