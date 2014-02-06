ifeq ($(DEBUG),)
LATEX = latex -interaction=nonstopmode 2>&1 >/dev/null
DVIPS = dvips -q 2>&1 >/dev/null
else
LATEX = latex -interaction=nonstopmode 
DVIPS = dvips
endif

OBJECTS=$(shell ls *.tex)
FILELIST=$(basename $(OBJECTS))

export TEXINPUTS := ..:
export DVIPSHEADERS := ..:

.PHONY: all clean veryclean

.SECONDARY: $(addsuffix .ps, $(FILELIST))

all: $(addsuffix .pdf, $(FILELIST))

%.pdf: %.ps
	ps2pdf $<

%.ps: %.dvi
	$(DVIPS) $<

%.dvi: %.tex
	$(LATEX) $<

clean: 
	$(RM) $(addsuffix .aux, $(FILELIST)) \
	      $(addsuffix .log, $(FILELIST)) \
	      $(addsuffix .dvi, $(FILELIST)) \
	      $(addsuffix .out, $(FILELIST))

veryclean: clean
	$(RM) $(addsuffix .ps, $(FILELIST)) \
	      $(addsuffix .pdf, $(FILELIST))
