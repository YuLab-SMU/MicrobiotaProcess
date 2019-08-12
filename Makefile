PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd pdf check clean

rd:
	Rscript -e 'library(methods); devtools::document()'

pdf: rd
	R CMD Rd2pdf --no-preview --force ../$(PKGSRC)
        
build: pdf 
	cd ..;\
    R CMD build $(PKGSRC)

check: build
	cd ..;\
    R CMD check --as-cran $(PKGNAME)_$(PKGVERS).tar.gz

install: check
	cd ..;\
    R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

clean:
	cd ..;\
    rm -rf $(PKGNAME).Rcheck
