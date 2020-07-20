ALL: ms.html

ms.html: ms.Rmd
	Rscript -e 'rmarkdown::render("$(<F)", output_format = "bookdown::html_document2")'
	open $@

%.pdf: %.Rmd
	Rscript -e 'rmarkdown::render("$(<F)", output_format = "bookdown::pdf_document2")'
	open $@
	
%.docx: %.Rmd
	Rscript -e 'rmarkdown::render("$(<F)", output_format = "bookdown::word_document2")'
	open $@

clean:
	rm ms.html ms.docx

diff.pdf: diff.tex 
	pdflatex -interaction nonstopmode diff
	# bibtex diff
	pdflatex -interaction nonstopmode diff
	open diff.pdf

	# diff.tex: geb_1st_submit.tex cooc.tex
	# 	latexdiff geb_1st_submit.tex cooc.tex > diff.tex

diff.tex: ms.tex ms_submitted.tex
	latexdiff ms_submitted.tex ms.tex > diff.tex
