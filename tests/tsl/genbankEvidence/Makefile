MODROOT = ../../../..

include ${MODROOT}/defs.mk

include ${ROOT}/make/rules.mk
include ${ROOT}/make/exprDefs.mk

gencodeOrg = human
ensemblVer = 71_37
ensemblCDnaDb = homo_sapiens_cdna_${ensemblVer}
grcRefAssembly = GRCh37
ucscGenomeDb = hg19

ensToUcscMappingPsl = output/ensToUcscMapping.psl

test:: mrnaEvidTest estEvidTest estAllEvidTest estInclSingleEvidTest ensemblChrMTest ensemblChr22Test ensemblIdListTest ensemblChr6HoxTest

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

ensemblChrMTest: mkdirs ${ensToUcscMappingPsl}
	${cdnaGetEnsemblAligns} -chrom=chrM -genomeDb=${ucscGenomeDb} ${ensemblCDnaOpts} ${ensemblCDnaDb} ${ensToUcscMappingPsl} output/$@.psl output/$@.missing.acc
	sort output/$@.psl > output/$@.sort.psl
	${cdnaPslEvidence} -ignoreMatch output/$@.sort.psl ${ucscGenomeSeqs} output/$@.alninter
	diff expected/$@.sort.psl output/$@.sort.psl
	diff expected/$@.alninter output/$@.alninter

ensemblChr22Test: mkdirs ${ensToUcscMappingPsl}
	${cdnaGetEnsemblAligns} -chrom=chr22 -genomeDb=${ucscGenomeDb} ${ensemblCDnaOpts} ${ensemblCDnaDb} ${ensToUcscMappingPsl} output/$@.psl output/$@.missing.acc

ensemblIdListTest: mkdirs ${ensToUcscMappingPsl}
	${cdnaGetEnsemblAligns} -accverList=input/ensembl-cdna.acc -genomeDb=${ucscGenomeDb} ${ensemblCDnaOpts} ${ensemblCDnaDb} ${ensToUcscMappingPsl} output/$@.psl output/$@.missing.acc
	sort output/$@.psl > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

# test mapping an MHC region
ensemblChr6HoxTest: mkdirs ${ensToUcscMappingPsl}
	${cdnaGetEnsemblAligns} -chrom=chr6_cox_hap2 -genomeDb=${ucscGenomeDb} ${ensemblCDnaOpts} ${ensemblCDnaDb} ${ensToUcscMappingPsl} output/$@.psl output/$@.missing.acc
	sort output/$@.psl > output/$@.sort.psl
	diff expected/$@.sort.psl output/$@.sort.psl

# mapping PSLs for cdnaGetEnsemblAligns
${ensToUcscMappingPsl}: mkdirs ${chrMPsl}
	${ensToUcscChromMap} ${ensemblCDnaDb} ${grcRefAssembly} ${ucscGenomeDb} $@.${tmpExt}
	mv -f $@.${tmpExt} $@

mkdirs:
	@mkdir -p output

clean::
	rm -rf output
