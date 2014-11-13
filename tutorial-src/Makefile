#
# Generates the phenomenological model tutorial
#
# The markdown document
MDFILE = phenomenological_model_tutorial.md

# all output formats
PDFDOC = $(MDFILE:.md=.pdf)
PDFSADOC = $(MDFILE:.md=_standalone.pdf)
IPYNBDOC = $(MDFILE:.md=.ipynb)
IPYNBBUILTDOC = $(MDFILE:.md=_built.ipynb)
HTMLBASICDOC = $(MDFILE:.md=_basic.html)
HTMLDOC = $(MDFILE:.md=.html)
WEBSITE_TARBALL = asymmetric-motion-website-files.tar.gz

MDPP = $(shell command -v ./mdpp || command -v mdpp)
MD2IPYNB = $(shell command -v ./md2ipynb || command -v md2ipynb)

PANDOC_VERSION = $(shell pandoc --version)

ifeq ($(MDPP),)
$(error mdpp not found: "wget https://raw.githubusercontent.com/ap--/mdpp/master/mdpp")
endif
ifeq ($(MD2IPYNB),)
$(error md2ipynb not found: "wget https://raw.githubusercontent.com/ap--/md2ipynb/master/md2ipynb")
endif


# IPYTHON NOTEBOOK OPTIONS
IPYMDPPARGS = --without .OPTIONALFIGURE --clean
PDFSAMDPPARGS = --without .INCLUDECODESNIPPET --clean
PDFMDPPARGS = --without .OPTIONALFIGURE --without .INCLUDECODESNIPPET --clean


# output format specific changes
PDFTEMPLATE = templates/current-biology-pdf-template.latex

# Get all relevant information from md file
PYCODEREGEXP = "s/<\!--@\s\+.INCLUDECODE\(SNIPPET\)\?\s\+python\s\+\(\S*\)\(\s\+\S*\)*@-->/\2/p" 

INCLUDEDPYFILES = $(shell sed -n $(PYCODEREGEXP) $(MDFILE) | sort -u)

SVGFIGURES = $(shell $(MDPP) $(MDFILE) \
	     | pandoc --default-image-extension=svg -f markdown -t json - \
	     | jq -r 'def recurse(f): ., (f | select(. != null) | recurse(f)); \
		      recurse(select(.>=[]) | .[]) | select(type|.=="object") |\
		      select(.t=="Image") | .c[1][0]')
COMMONFIGURES = $(shell $(MDPP) $(IPYMDPPARGS) $(MDFILE) \
	     | pandoc --default-image-extension=svg -f markdown -t json - \
	     | jq -r 'def recurse(f): ., (f | select(. != null) | recurse(f)); \
		      recurse(select(.>=[]) | .[]) | select(type|.=="object") |\
		      select(.t=="Image") | .c[1][0]')
INCLUDEDPYFIGURES = $(filter-out $(COMMONFIGURES),$(SVGFIGURES))


# get all figures from the markdown document
PDFFIGURES = $(SVGFIGURES:.svg=.pdf)

COMMONARGS = -f markdown+yaml_metadata_block -s

# PDF options
PDFARGS = -t latex \
	  --default-image-extension=pdf \
	  --template $(PDFTEMPLATE) \
	  --listings
# PDFVARS = -V documentclass=scrreprt
PDFVARS = -V papersize=a4paper \
	  -V fontsize=10pt


.PHONY: debug clean help

help:
	@echo "Supplemental Experimental Procedures Makefile:"
	@echo "----------------------------------------------"
	@echo "Generate the pdf as a standalone document:"
	@echo "> make pdf-standalone"
	@echo ""
	@echo "Generate the ipython notebook:"
	@echo "> make ipynb"
	@echo ""
	@echo "Generate a static, includable .html version of the ipython notebook:"
	@echo "> make html-basic"
	@echo ""
	@echo "Generate a full .html version of the ipython notebook:"
	@echo "> make html"
	@echo ""



all: pdf-standalone ipynb html-basic html

pdf-standalone: $(PDFSADOC)

ipynb: $(IPYNBDOC)

html-basic: $(HTMLBASICDOC)

html: $(HTMLDOC)

website-tarball: $(WEBSITE_TARBALL)

$(WEBSITE_TARBALL): $(IPYNBBUILTDOC) $(HTMLBASICDOC) $(HTMLDOC)
	tar -czf $@ $(IPYNBBUILTDOC) $(HTMLBASICDOC) $(HTMLDOC)

.PRECIOUS: %.svg %.pdf

$(INCLUDEDPYFIGURES): $(INCLUDEDPYFILES)
	for PY in "$(INCLUDEDPYFILES)"; do python $$PY; done

%.svg: %.py
	python $<

%.pdf: %.svg
	inkscape $< --export-pdf=$@

$(IPYNBDOC): $(MDFILE) $(INCLUDEDPYFILES) $(COMMONFIGURES) Makefile
	$(MDPP) $(IPYMDPPARGS) $< | $(MD2IPYNB) --out $@ -

$(IPYNBBUILTDOC): $(IPYNBDOC)
	runipy $< $@

$(PDFDOC): $(MDFILE) $(PDFFIGURES) $(PDFTEMPLATE) $(INCLUDEDPYFILES) Makefile
	$(MDPP) $(PDFMDPPARGS) $< | pandoc $(PDFARGS) $(PDFVARS) $(COMMONARGS) -o $@ -

$(PDFSADOC): $(MDFILE) $(PDFFIGURES) $(PDFTEMPLATE) $(INCLUDEDPYFILES) Makefile
	$(MDPP) $(PDFSAMDPPARGS) $< | pandoc $(PDFARGS) $(PDFVARS) $(COMMONARGS) -o $@ -

$(HTMLBASICDOC): $(IPYNBBUILTDOC)
	ipython nbconvert --to html --template basic $< --stdout > $@

$(HTMLDOC): $(IPYNBBUILTDOC)
	ipython nbconvert --to html --template full $< --stdout > $@

debug:
	echo "II: pandoc should be >= 1.11"
	pandoc --version
	echo ""

clean:
	rm -f $(SVGFIGURES)
	rm -f $(PDFFIGURES)
	rm -f $(PDFDOC)
	rm -f $(PDFSADOC)
	rm -f $(IPYNBDOC)
