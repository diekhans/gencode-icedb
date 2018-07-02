include ${ROOT}/tests/testDefs.mk

hsGencodeVer = V28
hsEnsemblVer = 92_38
hsEnsemblCDnaDb = homo_sapiens_cdna_${hsEnsemblVer}
hsGrcRefAssembly = GRCh38
hsGrcRefAssemblyReport = ${ROOT}/tests/tsl/etc/GCF_000001405.36_GRCh38.p10_assembly_report.txt
hsUcscGenomeDb = hg38

mmGencodeVer = VM16
mmEnsemblVer = 91_38
mmEnsemblCDnaDb = mus_musculus_cdna_${hsEnsemblVer}
mmGrcRefAssembly = GRCh38
mmGrcRefAssemblyReport = ${ROOT}/tests/tsl/etc/GCF_000001405.36_GRCh38.p10_assembly_report.txt
mmUcscGenomeDb = mm10


gencodeDbLoad = ${BINDIR}/gencodeDbLoad
tslGetUcscRnaAligns = ${BINDIR}/tslGetUcscRnaAligns
tslGetEnsemblRnaAligns = ${BINDIR}/tslGetEnsemblRnaAligns
tslGbffGetProblemCases = ${BINDIR}/tslGbffGetProblemCases
tslGenbankProblemCasesLoad = ${BINDIR}/tslGenbankProblemCasesLoad
tslGencodeCollectSupport = ${BINDIR}/tslGencodeCollectSupport
tslGencodeCollectSupportMkJobs = ${BINDIR}/tslGencodeCollectSupportMkJobs
tslGencodeCollectSupportFinishJobs = ${BINDIR}/tslGencodeCollectSupportFinishJobs

# use public database if not on hgwdev
ifneq (${hostname}, hgwdev)
   export HGDB_CONF=${ROOT}/tests/tsl/etc/hg.pub.conf
endif

# table names (match py definition)
UCSC_RNA_ALN_TBL = ucsc_rna_aln
UCSC_EST_ALN_TBL = ucsc_est_aln
ENSEMBL_RNA_ALN_TBL = ensembl_rna_aln
