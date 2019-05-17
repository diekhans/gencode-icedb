ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/rsl/rslTestDefs.mk

gencodeTestDb = output/gencode.db
sjTestDb =  output/sj.db

dbinput = ../evidence/input
dbexpect = ../evidence/expected

test:	supportTests \
	novelTests

##
# Test of support classify
##
supportTests:	supportCollectGenesTest \
		supportCollectTransTest \
		supportCollectMkJobsTest

supportCollectGenesTest: ${gencodeTestDb} ${sjTestDb} mkdirs
	${rslGencodeCollectSupport} ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.supp.tsv ENSG00000235169.7
	${diff} expected/$@.supp.tsv output/$@.supp.tsv

supportCollectTransTest: ${gencodeTestDb} ${sjTestDb} mkdirs
	${rslGencodeCollectSupport} ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.supp.tsv ENST00000444870.6 ENST00000462356.5
	${diff} expected/$@.supp.tsv output/$@.supp.tsv

supportCollectMkJobsTest: mkdirs
	rm -rf output/$@.work
	${rslGencodeCollectSupportMkJobs} --maxgenes=10 ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.work
	sed -Ee 's#^.*(/bin/rslGencodeCollectSupport)#??\1#' output/$@.work/batch.jobs > output/$@.batch.jobs.cleaned
	${diff} expected/$@.batch.jobs.cleaned output/$@.batch.jobs.cleaned
	${jobsToMake} output/$@.work/batch.jobs output/$@.mk
	${MAKE} -f output/$@.mk
	ls -1 output/$@.work/results/ > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${rslGencodeCollectSupportFinishJobs} output/$@.work output/$@.db
	$(call sqldumpdiff,gencode_intron_support)

##
# test of novel
##
novelTests: novelRange1Test \
	novelMkJobsTest

novelRange1Test: mkdirs
	${rslGencodeCollectNovel} ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.novel.tsv chr1:2193910-7731544
	${diff} expected/$@.novel.tsv output/$@.novel.tsv

# make sure putative introns start in window
novelSplitTest: mkdirs
	${rslGencodeCollectNovel} ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.win1.novel.tsv chr22:0-10940708
	${rslGencodeCollectNovel} ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.win2.novel.tsv chr22:10940708-20000000
	${diff} expected/$@.win1.novel.tsv output/$@.win1.novel.tsv
	${diff} expected/$@.win2.novel.tsv output/$@.win2.novel.tsv

novelMkJobsTest: mkdirs
	rm -rf output/$@.work
	${rslGencodeCollectNovelMkJobs} --testRange=chr22:10640100-17601091 --windowSize=1000000 ${genomeSeqsSpecs} ${gencodeTestDb} ${sjTestDb} output/$@.work
	sed -Ee 's#^.*(/bin/rslGencodeCollectNovel)#??\1#' output/$@.work/batch.jobs > output/$@.batch.jobs.cleaned
	${diff} expected/$@.batch.jobs.cleaned output/$@.batch.jobs.cleaned
	${jobsToMake} output/$@.work/batch.jobs output/$@.mk
	${MAKE} -f output/$@.mk
	ls -1 output/$@.work/results/ > output/$@.results.lst
	${diff} expected/$@.results.lst output/$@.results.lst
	${rslGencodeCollectNovelFinishJobs} output/$@.work output/$@.db
	$(call sqldumpdiff,gencode_intron_novel)

# create test input databases
${gencodeTestDb}:
	@mkdir -p $(dir $@)
	@rm -f $@.tmp
	${ucscGencodeDbLoad} $@.tmp --attrs=${dbinput}/set1.hs.attrs.v26.tsv --transcriptSource=${dbinput}/set1.hs.trans-source.v26.tsv --transcriptionSupportLevel=${dbinput}/set1.hs.trans-support-level.v26.tsv --tags=${dbinput}/set1.hs.trans-tag.v26.tsv --genes=${dbinput}/set1.hs.comp.v26.gp --pseudoGenes=${dbinput}/set1.hs.pseudo.v26.gp
	mv -f $@.tmp $@

${sjTestDb}:
	@mkdir -p $(dir $@)
	@rm -f $@.tmp
	${sraRunInfoDbLoad} $@.tmp ${dbinput}/set1.hs.run-info.tsv
	${rslMappingMetadataDbLoad} $@.tmp ${dbinput}/set1.hs.mapping-params.tsv ${dbinput}/set1.hs.stardir.tsv
	bgzip -c ${dbexpect}/mkStarMergeTest.merged.sjsup > $(basename $@).sjsup.gz
	tabix -0 -p bed $(basename $@).sjsup.gz
	mv -f $@.tmp $@


mkdirs:
	@mkdir -p output

clean:
	rm -rf output
