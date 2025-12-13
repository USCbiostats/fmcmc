
install:
	Rscript -e "devtools::install()"
		
build:
	R CMD build .

README.md: README.qmd
	quarto render README.qmd

.PHONY: checfull checkv clean

check:
	Rscript -e "devtools::check()"

man: R/* 
	Rscript --vanilla -e 'devtools::document()'

.PHONY: man

