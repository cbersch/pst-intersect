.SUFFIXES : .tex .ltx .dvi .ps .pdf .eps

PACKAGE = pst-intersect

LATEX = latex

ARCHNAME = $(PACKAGE)-$(shell date +"%y%m%d")
ARCHNAME_TDS = $(PACKAGE).tds

ARCHFILES = $(PACKAGE).dtx $(PACKAGE).ins Makefile \
            README Changes $(PACKAGE).pdf

PS2PDF = GS_OPTIONS=-dPDFSETTINGS=/prepress ps2pdf

all : doc

doc : $(PACKAGE).pdf

dist : doc Changes
	mkdir -p $(PACKAGE)
	cp $(ARCHFILES) $(PACKAGE)
	tar chvzf $(ARCHNAME).tar.gz $(PACKAGE)
	rm -rf $(PACKAGE)
	@ echo
	@ echo $(ARCHNAME).tar.gz

$(PACKAGE).dvi: L = english
$(PACKAGE)-DE.dvi: L = ngerman
%.dvi: $(PACKAGE).dtx $(PACKAGE).sty $(PACKAGE).ist $(PACKAGE).pro
	$(LATEX) -jobname=$(basename $@) '\newcommand*{\mainlang}{$(L)}\input{$(PACKAGE).dtx}'
	$(LATEX) -jobname=$(basename $@) '\newcommand*{\mainlang}{$(L)}\input{$(PACKAGE).dtx}'

%.ps: %.dvi
	dvips $< 
%.pdf: %.ps
	$(PS2PDF) $< $@
	python $<

$(PACKAGE).sty $(PACKAGE).pro $(PACKAGE).tex: $(PACKAGE).ins $(PACKAGE).dtx
	tex $<

Changes: Changes.py $(PACKAGE).dtx
	python $<

arch : Changes
	zip $(ARCHNAME).zip $(ARCHFILES)

arch-tds : Changes
	$(RM) $(ARCHNAME_TDS).zip
	mkdir -p tds/tex/latex/$(PACKAGE)
	mkdir -p tds/tex/generic/$(PACKAGE)
	mkdir -p tds/doc/latex/$(PACKAGE)
	mkdir -p tds/source/latex/$(PACKAGE)
	mkdir -p tds/dvips/$(PACKAGE)
	cp $(PACKAGE).sty tds/tex/latex/$(PACKAGE)/
	cp $(PACKAGE).tex tds/tex/generic/$(PACKAGE)/
	cp $(PACKAGE).pro tds/dvips/$(PACKAGE)/
	cp Changes $(PACKAGE).pdf README tds/doc/latex/$(PACKAGE)/
	cp $(PACKAGE).dtx $(PACKAGE).ins Makefile \
	  tds/source/latex/$(PACKAGE)/
	cd tds ; zip -r ../$(ARCHNAME_TDS).zip tex doc source dvips
	cd ..
	rm -rf tds

ctan : dist arch-tds
	tar cf $(PACKAGE).tar $(ARCHNAME_TDS).zip $(ARCHNAME).tar.gz

clean :
	$(RM) $(addprefix $(PACKAGE), .dvi .ps .log .aux .bbl .blg .out .tmp .toc .hd))

veryclean : clean
	$(RM) $(addprefix $(PACKAGE), .pdf .tex .sty .pro) Changes
