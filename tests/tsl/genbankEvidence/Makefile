ROOT = ../../..

binDir = ${ROOT}/bin
getUcscRnaAligns = ${binDir}/getUcscRnaAligns
getEnsemblRnaAligns = ${binDir}/getEnsemblRnaAligns

gencodeOrg = human
ensemblVer = 87_38
ensemblCDnaDb = homo_sapiens_cdna_${ensemblVer}
grcRefAssembly = GRCh38
grcRefAssemblyReport = input/GCF_000001405.36_GRCh38.p10_assembly_report.txt
ucscGenomeDb = hg38

test:: ucscTests ensemblTests

##
# UCSC tsts
##
ucscTests: ucscRnaChrMTest ucscEstChrMTest

ucscRnaChrMTest: mkdirs
	rm -f output/$@.db
	${getUcscRnaAligns} -chrom=chrM ${ucscGenomeDb} rna output/$@.db ucscRnaAligns
	sqlite3 --batch output/$@.db 'select * from ucscRnaAligns;' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

ucscEstChrMTest: mkdirs
	rm -f output/$@.db
	${getUcscRnaAligns} -chrom=chrM ${ucscGenomeDb} est output/$@.db ucscEstAligns
	sqlite3 --batch output/$@.db 'select * from ucscEstAligns;' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

##
# Ensembl
##n
ensemblTests: ensemblChrMTest ensemblIdListTest

ensemblChrMTest: mkdirs
	@rm -f output/$@.db
	${getEnsemblRnaAligns} --ensChrom=MT ${ensemblCDnaDb} ${grcRefAssemblyReport} output/$@.db ensemblCDnaAln
	sqlite3 --batch output/$@.db 'select * from ensemblCDnaAln;' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

ensemblIdListTest: mkdirs
	@rm -f output/$@.db
	${getEnsemblRnaAligns} --accverList=input/ensembl-cdna.acc ${ensemblCDnaDb} ${grcRefAssemblyReport} output/$@.db ensemblCDnaAln
	sqlite3 --batch output/$@.db 'select * from ensemblCDnaAln;' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
