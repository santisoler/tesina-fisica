# Build the PDF of the manuscript.
#
# Commands:
#
#	make            Compile PDF and cleans the auxiliary files
#	make build		Only compiles the PDF
#	make clean  	Delete all auxiliary files
#   make clean_pdf  Delete PDF file

##############################################################################

# Configuration

PROJECT = tesina-soler
SRC = .
BIB = bibtex
STYLES = bibtex

## Latex, bibtex and styles files
TEX_FILES = $(shell find $(SRC) -maxdepth 1 -name '*.tex')
BIB_FILES = $(shell find $(BIB) -maxdepth 1 -name '*.bib')
BST_FILES = $(shell find $(STYLES) -maxdepth 1 -name '*.bst')

### Compilation Flags
LATEX_FLAGS  = -halt-on-error


##############################################################################

# Targets

all: build clean

build: $(SRC)/$(PROJECT).pdf

clean_pdf:
	@ echo "Deleting PDF"
	@ rm $(SRC)/$(PROJECT).pdf

clean:
	@ echo "Cleaning aux files..."
	@ rm *.aux
	@ rm */*.aux
	@ rm *.toc
	@ rm *.bbl
	@ rm *.blg
	@ rm *.log

$(SRC)/$(PROJECT).aux: $(TEX_FILES) $(STY_FILES)
	pdflatex $(LATEX_FLAGS) $(SRC)/$(PROJECT) >/dev/null

$(SRC)/$(PROJECT).bbl: $(BIB_FILES) $(BST_FILES) | $(SRC)/$(PROJECT).aux
	bibtex $(PROJECT) >/dev/null
	pdflatex $(LATEX_FLAGS) $(SRC)/$(PROJECT) >/dev/null

$(SRC)/$(PROJECT).pdf: $(SRC)/$(PROJECT).aux $(SRC)/$(PROJECT).bbl
	pdflatex $(LATEX_FLAGS) $(SRC)/$(PROJECT) >/dev/null
	
