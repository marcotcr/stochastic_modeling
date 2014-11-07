all:
	R -e "library(knitr); knit('hw.Rnw')"
	pdflatex hw.tex

clean:
	rm -r *.aux *.log figure cache hw.tex hw.pdf
