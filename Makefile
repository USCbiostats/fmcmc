
install:
	Rscript -e "devtools::install()"
		
build:
	R CMD build .

inst/NEWS: NEWS.md
	Rscript -e "rmarkdown::pandoc_convert('NEWS.md', 'plain', output='inst/NEWS')"&& \
	head -n 80 inst/NEWS

README.md: README.qmd
	quarto render README.qmd

.PHONY: checfull checkv clean

check:
	Rscript -e "devtools::check()"

man: R/* 
	Rscript --vanilla -e 'devtools::document()'

.PHONY: man

