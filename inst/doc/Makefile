


hlookup.ps: hlookup.tex
	latex hlookup; latex hlookup; dvips -o hlookup.ps hlookup

hlookup.tex: hlookup.Rnw
	sweave hlookup.Rnw

clean:
	rm -f hlookup.tex hlookup-*.eps hlookup.dvi hlookup.log hlookup.aux
