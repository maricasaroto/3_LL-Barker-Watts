WOS := $(shell uname)

DOC := $(shell egrep -l '^[^%]*\\begin\{document\}' *.tex | cut -f1 -d.)

BIB := $(shell egrep -l '^[^%]*\\bibliography\{' $(DOC).tex)

pdflatex:
	@pdflatex -shell-escape $(DOC).tex
ifeq ($(BIB), )
	@pdflatex -shell-escape $(DOC).tex
else
	@bibtex $(DOC)
	@pdflatex -shell-escape $(DOC).tex
	@pdflatex -shell-escape $(DOC).tex
endif
ifeq ($(WOS),Darwin)
	@open $(DOC).pdf
else
	@evince $(DOC).pdf &
endif

clean:
	@rm -f *.aux *.bbl *.blg *.log *.out *.toc *.dvi *.nav *.snm *.pyg
	@rm -f tex/*.aux 
