.PHONY: Main.pdf pdf clean distclean

Main.pdf: Main.tex Bibliografia/Bib.bib 
		@pdflatex Main
		@bibtex Main
		@pdflatex Main-frn
		@pdflatex  Main
#		texi2pdf main

pdf: Main.pdf

clean: 
	rm -f *.aux *.log *.out
        
distclean: 
	rm -f *.aux *.log *.out *.pdf *.gz
