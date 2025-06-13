all: README.pdf mini-papers.pdf references/README.md

clean: mini-papers-clean graphs-clean

README.pdf: graphs-all

################################################################################
# Render graphs and re-create data if necessary

GRAPH_SCRIPTS = $(wildcard graphs/*.py)
.PHONY: graphs-all
graphs-all: $(GRAPH_SCRIPTS:.py=.svg)

graphs/%.png: graphs/%.py
	python3 $<

graphs/%.svg: graphs/%.py
	python3 $<

graphs/%.small.svg: graphs/%.svg
	scour $< $@ --shorten-ids --no-line-breaks --enable-comment-stripping --remove-metadata --enable-id-stripping

graphs-clean:
	rm -f graphs/*.png graphs/*.svg

################################################################################
# General typesetting & references

%.pandoc.pdf: %.md typeset/header.tex typeset/ieee.csl references/ref.bib Makefile
	pandoc --bibliography=references/ref.bib --csl=typeset/ieee.csl -L typeset/columns.lua --metadata numbersections=true  --metadata link-citations=true --metadata ragged-columns=true -o "$@" --include-in-header=typeset/header.tex "$<"

%.linear.pdf: %.pandoc.pdf
	qpdf --linearize --remove-unreferenced-resources=yes $< $@

%.pdf: %.linear.pdf
	gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dSAFER -dNODOCINFO -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH -sOutputFile=$@ $<

MINIPAPER_PDFS := $(wildcard mini-papers/*.md)

mini-papers-clean:
	rm -f mini-papers/*.pdf

.PHONY: mini-papers-all
mini-papers.pdf: $(MINIPAPER_PDFS:.md=.pdf)
	qpdf --empty --pages $^ -- $@

# Verifies all local files that are referenced exist in our references dir
check-ref-files:
	cat references/ref.bib | grep file | grep -v 'url =' | cut -f 2 -d '{' | cut -f 1 -d '}' | xargs -n 1 echo | awk -e '{print "references/" $$1 }' | xargs ls -lh

references/README.md: references/ref.bib
	cd references && python3 generate_readme.py

################################################################################
# Re-create the data, for given bit-sizes, also compress the dataset

data/%.sqlite3.xz:
	xz -k -vvz -T 0 "data/$*.sqlite3"

.PHONY: data/%.sqlite3
data/%.sqlite3:
	python3 steps/1-primes.py $*
	python3 steps/2-cornacchia.py $*
	python3 steps/3-trial-division.py $*
	python3 steps/4-generator.py $*
	python3 steps/5-curves.py $*
	python3 steps/6-glv.py $*
	python3 steps/7-curvefactor.py $*
	python3 steps/8-embedding.py $*
