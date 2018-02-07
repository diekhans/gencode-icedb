##
# Semi-manual testing that builds a full database.  Used in development
# build not for regression, as it is slow.
##
ROOT = ../../..
include ${ROOT}/config.mk
include ${ROOT}/tests/tsl/tslTestDefs.mk

evidDir = evid.d
evidDb = ${evidDir}/hsEvid.db

annonDir = annon.d
annonDb = ${annonDir}/hsAnnon.db
inDbs = ${evidDb} ${annonDb}


test::



##
# Building of evidence db (~1.25 hours)
##

${evidDb}:
	@mkdir -p $(dir $@)
	rm -f ${evidDb} ${evidDb}.tmp
	${tslGetEnsemblRnaAligns} ${hsEnsemblCDnaDb} ${hsGrcRefAssemblyReport} ${evidDb}.tmp ${ENSEMBL_RNA_ALN_TBL}
	${tslGetUcscRnaAligns} ${hsUcscGenomeDb} rna ${evidDb}.tmp ${UCSC_RNA_ALN_TBL}
	${tslGetUcscRnaAligns} ${hsUcscGenomeDb} est ${evidDb}.tmp ${UCSC_EST_ALN_TBL}
	mv -f ${evidDb}.tmp ${evidDb}
##
# building annotation db (~20 seconds)
##
${annonDb}:
	@mkdir -p $(dir $@)
	rm -f ${annonDb} ${annonDb}.tmp
	${gencodeDbLoad} --hgdb=${hsUcscGenomeDb} --version=${hsGencodeVer} ${annonDb}.tmp
	mv -f ${annonDb}.tmp ${annonDb}

clean::
	rm -rf output

realclean: clean
	rm -rf ${evidDir} ${annonDir}
