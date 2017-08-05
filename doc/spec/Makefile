ROOT = .
.SECONDARY:

outdir = output

docs += gencode-icedb-spec

docsPdf = ${docs:%=${outdir}/%.pdf}
docsFontchk = ${docs:%=${outdir}/%.fontchk}
drawings = 
#${outdir}/ice-db-arch.pdf
images = $(wildcard images/*.pdf)
depends = ${drawings} ${images}

all:: ${docsPdf} ${docsFontchk} ${depends}

ltxopts = -file-line-error-style -output-directory=output
ltxmode =\\nonstopmode\\input

${outdir}/%.pdf: %.ltx ${depends}
	@mkdir -p $(dir $@)
	pdflatex ${ltxopts} -jobname $*.tmp ${ltxmode} $< </dev/null || (rm -f $@; false)
	bibtex ${outdir}/$*.tmp </dev/null || (rm -f $@; false)
	pdflatex ${ltxopts} -jobname $*.tmp ${ltxmode} $< </dev/null || (rm -f $@; false)
	pdflatex ${ltxopts} -jobname $*.tmp ${ltxmode} $< </dev/null || (rm -f $@; false)
	mv -f ${outdir}/$*.tmp.pdf $@

${outdir}/%.pdf: drawings/%.svg
	@mkdir -p $(dir $@)
	inkscape --export-pdf=$@.tmp $<
	mv -f $@.tmp $@

${outdir}/%.fontchk: ${outdir}/%.pdf
	pdfFontChk $<
	touch $@

${outdir}/%.pdf: drawings/%.svg
	@mkdir -p $(dir $@)
	inkscape --export-pdf=$@.tmp $<
	mv -f $@.tmp $@

${outdir}/%.pdf: ${outdir}/%.eps
	epspdf $< $@.tmp.pdf
	mv -f $@.tmp.pdf $@

${outdir}/%.eps: drawings/%.gnuplot
	@mkdir -p $(dir $@)
	gnuplot $< >$@.tmp
	mv -f $@.tmp $@

clean::
	rm -rf output
