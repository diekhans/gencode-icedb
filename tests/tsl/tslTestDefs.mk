# force sort to be consistent
export LC_ALL=C
.SECONDARY:


gencodeDbLoad = ${BINDIR}/gencodeDbLoad
tslGetUcscRnaAligns = ${BINDIR}/tslGetUcscRnaAligns
tslGetEnsemblRnaAligns = ${BINDIR}/tslGetEnsemblRnaAligns
tslGbffGetProblemCases = ${BINDIR}/tslGbffGetProblemCases
tslGenbankProblemCasesLoad = ${BINDIR}/tslGenbankProblemCasesLoad
tslGencodeCollectSupport = ${BINDIR}/tslGencodeCollectSupport
tslGencodeCollectSupportMkJobs = ${BINDIR}/tslGencodeCollectSupportMkJobs
tslGencodeCollectSupportFinishJobs = ${BINDIR}/tslGencodeCollectSupportFinishJobs

hsGencodeVer = V27
hsEnsemblVer = 91_38
hsEnsemblCDnaDb = homo_sapiens_cdna_${hsEnsemblVer}
hsGrcRefAssembly = GRCh38
hsGrcRefAssemblyReport = ${ROOT}/tests/tsl/etc/GCF_000001405.36_GRCh38.p10_assembly_report.txt
hsUcscGenomeDb = hg38

# always use public database, since hgwdev might be ahead
export HGDB_CONF=${ROOT}/tests/tsl/etc/hg.pub.conf

# table names (match py definition)
UCSC_RNA_ALN_TBL = ucsc_rna_aln
UCSC_EST_ALN_TBL = ucsc_est_aln
ENSEMBL_RNA_ALN_TBL = ensembl_rna_aln
