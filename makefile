default:
	@echo targets:  roxy install test

all: roxy install test

roxy:
	R -e "devtools::document()"

install:
	R CMD INSTALL .

test:
	(cd inst/unitTests; make)

