PROLOGUES=$(shell kpsewhich pstricks.pro pst-func.pro) ../../pst-intersect.pro
FILELIST = $(basename $(wildcard pst-intersect-testgraph*.ps))

.PHONY: all clean

all: $(addsuffix .pdf, $(FILELIST))

%.pdf: %.ps
	gs -q -P- -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sstdout=%stderr -sOutputFile=$@ -c .setpdfwrite -f $(PROLOGUES) $<

clean:
	$(RM) -f $(addsuffix .pdf, $(FILELIST))
