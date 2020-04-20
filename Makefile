.PHONY: all clean

all: big_small_paper.pdf

big_small_paper.pdf : figures/patch_cartoon.pdf

clean: 
	-rm *.aux *.log *.lof *.lot *.fff *.ttt *.out *.bbl *.blg

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	-bibtex $<

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.tiff : %.png
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

%.fonts.pdf : %.pdf
	# this will embed fonts properly
	inkscape $< --export-pdf=$@
