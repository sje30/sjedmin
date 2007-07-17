all: hlookup.pdf clean

hlookup.pdf: hlookup.tex
	pdflatex $^

hlookup.ps: hlookup.tex
	latex hlookup; latex hlookup; dvips -o hlookup.ps hlookup

hlookup.tex: hlookup.Rnw
	echo 'Sweave("$^")' | R --vanilla -q


# Does not delete hlookup.pdf, final requirement.
.PHONY:
clean:
	rm -f hlookup.tex hlookup-*.pdf  \
	hlookup-*.eps hlookup.dvi hlookup.log hlookup.aux

