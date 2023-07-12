PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)
BIOCVER := RELEASE_3_17

all: rd check clean

alldocs: rd readme

rd:
	Rscript -e 'library(methods); devtools::document()'

readme:
	Rscript -e 'rmarkdown::render("README.Rmd")'

codemetar:
	Rscript -e 'codemetar::write_codemeta()'

build:
	cd ..;\
	R CMD build $(PKGSRC)

build2:
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install:
	cd ..;\
	R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz

check: rd build
	cd ..;\
	Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

check2: rd build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

check3: rd build2
	cd ..;\
	R CMD check --ignore-vignettes $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

bignore:
	Rscript -e 'usethis::use_build_ignore(glob2rx("inst/figures/*.png"), escape = FALSE)'
	Rscript -e 'usethis::use_build_ignore(c("Makefile", "README.md", "README.Rmd", "CONDUCT.md", ".Rproj.user", ".Rproj"))'

gpcheck:
	Rscript -e 'goodpractice::gp()'

debug: rd build2 install

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

clean2:
	cd ..;\
	$(RM) $(PKGNAME)_$(PKGVERS).tar.gz

gitmaintain:
	git gc --auto;\
	git prune -v;\
	git fsck --full

rmoldrelease:
	git branch -D $(BIOCVER)

release:
	git checkout $(BIOCVER):\
	git fetch --all

update:
	git fetch --all;\
	git checkout devel;\
	git merge upstream/devel;\
	git merge origin/devel

push: update
	git push upstream devel;\
	git push origin devel

biocinit:
	git remote add upstream git@git.bioconductor.org:packages/$(PKGNAME).git;\
	git fetch --all
