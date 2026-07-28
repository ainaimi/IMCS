# Introduction to Monte Carlo Simulation (IMCS)

Course materials for *Introduction to Monte Carlo Simulation*, taught by
[Ashley I. Naimi](https://bit.ly/3GFu2o1). The course is an applied
introduction to Monte Carlo simulation in the quantitative sciences: the role,
scope, and purpose of simulation; when and why to use it; its challenges and
limitations; and how to design, deploy, and analyze simulation studies in R.
By the end, students can implement their own simulation to estimate bias, mean
squared error, confidence interval coverage, and other statistics for an
estimator of their choice.

## Repository layout

```
lectures/                      the canonical course, one folder per section
  00_CourseIntro/                welcome and course logistics
  01_Section1_BackgroundContext/ why simulate?
  02_Section2_BasicTools/        functions; distributions; regression and distributions
  03_Section3_Design/            simulation design; Monte Carlo integration (ATT, CDE)
  04_Section4_PerformanceMeasures/ bias, MSE, coverage, and friends
  05_Section5_PresentingResults/ nested loop plots, zip plots
  06_Section6_SimulationInPractice/ worked studies: observational, RCT, mediation, CDM
  07_Section7_Computing/         parallel computing
  08_Section8_ComputingCluster/  cluster computing and SLURM
  appendix/                      simulation reading list
  misc/                          shared preamble.tex, style.css, and ref.bib used by every lecture
_images/                       figures used by the lectures (aRtsy.R generates the art)
syllabus.md                    course syllabus (Emory semester version)
install.R                      one-shot installer for the R packages the lectures use
```

The course is reformulated for different venues (semester course, SER and
Statistical Horizons workshops); this repository holds the single canonical
version. Delivered per-offering materials are archived outside the repo.

## Building a lecture

Requirements: R, a TeX distribution with `xelatex` (TeX Live, or
`tinytex::install_tinytex()`), and the R packages in `install.R` — each
lecture also auto-installs its own via `pacman::p_load()`.

Knit from inside the lecture's folder, e.g.:

```bash
cd lectures/04_Section4_PerformanceMeasures
Rscript -e "rmarkdown::render('section4_PerformanceMeasures.Rmd')"
```

Lectures use the tint (Tufte-style) PDF format. The shared preamble is
referenced relative to the lecture folder (`../misc/preamble.tex`), and figure
paths resolve from the repo root via `here()` (anchored by `IMCS.Rproj`).
Knitting a lecture also regenerates its companion `.R` script automatically (a
knitr purl hook), so each `.R` file is an extract of its `.Rmd`.

Long-running simulations are cached: section 6 reads its
`simulation_results*.csv` files at knit time instead of re-running the
simulations (the generating code is in the Rmd chunks).

## A note on PDFs

This repository tracks **source files only** — no compiled PDFs. Knitting
regenerates everything, with one caveat: five figure PDFs in `_images/` are
knit-time *inputs* that the lectures cannot regenerate themselves:

| figure | used by | how to obtain |
|---|---|---|
| `mediation_dag.pdf`, `rct_dag.pdf` | sections 3, 6, 8 | compile the `.tex` sources in `_images/` (xelatex, then pdfcrop) |
| `triangle_dag.pdf` | sections 3, 6 | no source in repo — contact the author |
| `par_ser.pdf` | section 7 | no source in repo — contact the author |
| `pareto_comparison.pdf` | section 2 (distributions) | no source in repo — contact the author |

Six other figure PDFs (`nested_loop_plot*`, `zip_plot_*`,
`weibull_exponential_comparison`) are written by the lectures themselves
during knit, so they appear after the first build.

## License

MIT — see `LICENSE.txt`.
