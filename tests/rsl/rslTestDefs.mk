include ${ROOT}/tests/testDefs.mk

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
