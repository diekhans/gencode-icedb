ROOT = ../../..

binDir = ${ROOT}/bin
cdnaEnsemblAligns = ${binDir}/cdnaEnsemblAligns2

gencodeOrg = human
ensemblVer = 87_38
ensemblCDnaDb = homo_sapiens_cdna_${ensemblVer}
grcRefAssembly = GRCh38
grcRefAssemblyReport = input/GCF_000001405.36_GRCh38.p10_assembly_report.txt
ucscGenomeDb = hg38

test:: mrnaEvidTest estEvidTest estAllEvidTest estInclSingleEvidTest

mrnaEvidTest: mkdirs
	${cdnaPslEvidence} input/mrnaEvid.psl ${ucscGenomeSeqs} output/$@.alninter
	diff expected/$@.alninter output/$@.alninter

estEvidTest: mkdirs
	${estPslEvidence} input/estEvid.psl ${ucscGenomeSeqs} output/$@.alninter
	diff expected/$@.alninter output/$@.alninter

estAllEvidTest: mkdirs
	${estPslEvidence} -allMultiExon input/estEvid.psl ${ucscGenomeSeqs} output/$@.alninter
	diff expected/$@.alninter output/$@.alninter

estInclSingleEvidTest: mkdirs
	${estPslEvidence} -inclSingleExon input/estEvid.psl ${ucscGenomeSeqs} output/$@.alninter
	diff expected/$@.alninter output/$@.alninter

##
# Ensembl
##
ensemblTests: ensemblChrMTest ensemblIdListTest

ensemblChrMTest: mkdirs
	@rm -f output/$@.db
	${cdnaEnsemblAligns} --ensChrom=MT ${ensemblCDnaDb} ${grcRefAssemblyReport} output/$@.db ensemblCDnaAln
	sqlite3 --batch output/$@.db 'select * from ensemblCDnaAln;' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

ensemblIdListTest: mkdirs
	@rm -f output/$@.db
	${cdnaEnsemblAligns} --accverList=input/ensembl-cdna.acc ${ensemblCDnaDb} ${grcRefAssemblyReport} output/$@.db ensemblCDnaAln
	sqlite3 --batch output/$@.db 'select * from ensemblCDnaAln;' | cut -f 2- | sort > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
