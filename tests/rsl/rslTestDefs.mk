export LC_ALL=C
.SECONDARY:

rslEncodeDccQuery = ${BINDIR}/rslEncodeDccQuery

rslSraRunInfoFilter = ${BINDIR}/rslSraRunInfoFilter
rslSraRunInfoDbLoad = ${BINDIR}/rslSraRunInfoDbLoad
rslMappingMetadataDbLoad = ${BINDIR}/rslMappingMetadataDbLoad

rslStarSjOutSplit = ${BINDIR}/rslStarSjOutSplit
rslMkStarSjOutSplits = ${BINDIR}/rslMkStarSjOutSplits
rslMkStarSjSupMergeJobs = ${BINDIR}/rslMkStarSjSupMergeJobs

gencodeDbLoad = ${BINDIR}/gencodeDbLoad
rslGencodeCollectSupport = ${BINDIR}/rslGencodeCollectSupport
rslGencodeCollectSupportMkJobs = ${BINDIR}/rslGencodeCollectSupportMkJobs
rslGencodeCollectSupportFinishJobs = ${BINDIR}/rslGencodeCollectSupportFinishJobs

rslGencodeCollectNovel = ${BINDIR}/rslGencodeCollectNovel
rslGencodeCollectNovelMkJobs = ${BINDIR}/rslGencodeCollectNovelMkJobs
rslGencodeCollectNovelFinishJobs = ${BINDIR}/rslGencodeCollectNovelFinishJobs

genomeSeqsHg38 = /hive/data/genomes/hg38/hg38.2bit
genomeSeqsSpecs = --genomeSeqs=${genomeSeqsHg38}

# remove {make out exists ... } from job file so it can be executed as commands
jobsToCmds = sed -Ee 's/\{check out exists (.+)\}/\1/'

# basic command to dump a database created by a test
sqldumpcmd = sqlite3 -header -batch output/$@.db

# call function with table as argument to dump table and diff with expected
sqldumpdiff = ${sqldumpcmd} 'select * from $(1)' > output/$@.$(1).tsv && ${diff} expected/$@.$(1).tsv output/$@.$(1).tsv

