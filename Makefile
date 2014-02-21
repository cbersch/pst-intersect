.SUFFIXES : .tex .ltx .dvi .ps .pdf .eps

PACKAGE = pst-intersect

LATEX = latex

ARCHNAME = $(PACKAGE)-$(shell date +"%y%m%d")
ARCHNAME_TDS = $(PACKAGE).tds

ARCHFILES = $(PACKAGE).dtx $(PACKAGE).ins Makefile \
            README Changes $(PACKAGE).pdf $(PACKAGE)-DE.pdf

PS2PDF = GS_OPTIONS=-dPDFSETTINGS=/prepress ps2pdf

all : doc doc-DE

doc : $(PACKAGE).pdf
doc-DE : $(PACKAGE)-DE.pdf

dist : doc doc-DE Changes
	mkdir -p $(PACKAGE)
	cp $(ARCHFILES) $(PACKAGE)

$(PACKAGE).dvi: L = english
$(PACKAGE)-DE.dvi: L = ngerman
%.dvi: $(PACKAGE).dtx $(PACKAGE).sty $(PACKAGE).tex $(PACKAGE).pro
	$(LATEX) -jobname=$(basename $@) '\newcommand*{\mainlang}{$(L)}\input{$(PACKAGE).dtx}'
	$(LATEX) -jobname=$(basename $@) '\newcommand*{\mainlang}{$(L)}\input{$(PACKAGE).dtx}'

%.ps: %.dvi
	dvips $< 
%.pdf: %.ps
	$(PS2PDF) $< $@

$(PACKAGE).sty $(PACKAGE).pro $(PACKAGE).tex: $(PACKAGE).ins $(PACKAGE).dtx
	tex $<

Changes: Changes.py $(PACKAGE).dtx
	python $<

arch-tds : Changes doc doc-DE
	$(RM) $(ARCHNAME_TDS).zip
	mkdir -p tds/tex/latex/$(PACKAGE)
	mkdir -p tds/tex/generic/$(PACKAGE)
	mkdir -p tds/doc/latex/$(PACKAGE)
	mkdir -p tds/source/latex/$(PACKAGE)
	mkdir -p tds/dvips/$(PACKAGE)
	cp $(PACKAGE).sty tds/tex/latex/$(PACKAGE)/
	cp $(PACKAGE).tex tds/tex/generic/$(PACKAGE)/
	cp $(PACKAGE).pro tds/dvips/$(PACKAGE)/
	cp Changes $(PACKAGE).pdf $(PACKAGE)-DE.pdf README tds/doc/latex/$(PACKAGE)/
	cp $(PACKAGE).dtx $(PACKAGE).ins Makefile \
	  tds/source/latex/$(PACKAGE)/
	cd tds ; zip -r ../$(ARCHNAME_TDS).zip tex doc source dvips
	cd ..
	rm -rf tds

ctan : dist arch-tds
	zip -r $(PACKAGE).zip $(ARCHNAME_TDS).zip $(PACKAGE)
	$(RM) -rf $(PACKAGE)/

clean :
	$(RM) $(foreach prefix, $(PACKAGE) $(PACKAGE)-DE, \
	        $(addprefix $(prefix), .dvi .ps .log .aux .bbl .blg .out .tmp \
	           .toc .hd))

veryclean : clean
	$(RM) $(addprefix $(PACKAGE), .pdf .tex .sty .pro) $(PACKAGE)-DE.pdf Changes
