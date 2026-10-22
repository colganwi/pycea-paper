# Uses tectonic from the pycea-paper conda env
main.pdf: main.tex sections/*.tex references.bib $(wildcard figures/*)
	PATH=/lab/solexa_weissman/wcolgan/tools/miniconda3/envs/pycea-paper/bin:$$PATH tectonic main.tex
clean:
	rm -f main.pdf *.aux *.bbl *.blg *.log
