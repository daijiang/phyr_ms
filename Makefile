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
