# Makefile for the IMCS course.
#
#   make            build every lecture PDF (and the syllabus)
#   make -k         same, but keep going past a failing lecture (full regression)
#   make <path.pdf> build a single lecture, e.g.:
#                     make lectures/04_Section4_PerformanceMeasures/section4_PerformanceMeasures.pdf
#   make figures    rebuild the DAG figures from their tex sources in _images/
#   make deps       install the R packages the lectures use (runs install.R)
#   make clean      remove LaTeX debris (aux/log/toc/synctex)
#   make clean-pdfs remove all generated PDFs (never the static figure inputs)
#
# Builds run serially. The old reason `make -j` was unsafe -- two section 6
# lectures racing to write the shared _images/nested_loop_plot_rct.pdf -- is
# gone now that only the RCT lecture writes it, but parallel builds remain
# untested.
# Knitting a lecture also regenerates its companion .R script (purl hook).

RSCRIPT := Rscript

LECTURE_RMDS := $(wildcard lectures/*/*.Rmd)
LECTURE_PDFS := $(LECTURE_RMDS:.Rmd=.pdf)

# figures the lectures regenerate themselves at knit time (cleanable)
GENERATED_FIGS := _images/nested_loop_plot.pdf _images/nested_loop_plot2.pdf \
                  _images/nested_loop_plot_rct.pdf _images/zip_plot_version1.pdf \
                  _images/zip_plot_version2.pdf _images/weibull_exponential_comparison.pdf

SHARED := lectures/misc/preamble.tex lectures/misc/ref.bib

.DELETE_ON_ERROR:
.PHONY: all figures deps clean clean-pdfs

all: $(LECTURE_PDFS) syllabus.pdf

# every lecture rebuilds when its source or the shared preamble/bibliography changes
lectures/%.pdf: lectures/%.Rmd $(SHARED)
	cd $(dir $@) && $(RSCRIPT) -e "rmarkdown::render('$(notdir $<)', quiet = TRUE)"

syllabus.pdf: syllabus.md
	$(RSCRIPT) -e "rmarkdown::render('syllabus.md', quiet = TRUE)"

# static figures that specific lectures include at knit time
lectures/02_Section2_BasicTools/section2_distributions.pdf: _images/pareto_comparison.pdf
lectures/03_Section3_Design/section3_SimulationDesign.pdf: _images/triangle_dag.pdf _images/mediation_dag.pdf
lectures/06_Section6_SimulationInPractice/section6_simulation.pdf: _images/triangle_dag.pdf
lectures/06_Section6_SimulationInPractice/section6_simulation_RCT.pdf: _images/rct_dag.pdf
lectures/06_Section6_SimulationInPractice/section6_simulation_mediation.pdf: _images/mediation_dag.pdf

# DAG figures compile from their tex sources
figures: _images/mediation_dag.pdf _images/rct_dag.pdf _images/triangle_dag.pdf

_images/%.pdf: _images/%.tex
	cd _images && xelatex -interaction=batchmode $(notdir $<) && rm -f $(notdir $(basename $<)).aux $(notdir $(basename $<)).log

# these ship with the course but have no generating source in the repo (see README)
_images/pareto_comparison.pdf:
	$(error $@ is a static figure input with no generating source in this repo; see the README)

deps:
	$(RSCRIPT) install.R

clean:
	find lectures _images \( -name "*.log" -o -name "*.aux" -o -name "*.toc" -o -name "*.synctex.gz" -o -name "tmp-pdfcrop-*.tex" \) -delete

clean-pdfs: clean
	rm -f $(LECTURE_PDFS) syllabus.pdf $(GENERATED_FIGS)
