# 🧬 BF591: R for Biological Sciences

> Structured, reproducible R workflows for biological data analysis

> Boston University – Bioinformatics Program

---

## Overview

This repository contains a structured series of assignments from **BF591: R for Biological Sciences**, focused on building reproducible, testable, and modular workflows for biological data analysis in R.

The course emphasizes:

* Functional programming in R
* Reproducible research practices
* Data wrangling with tidyverse
* Statistical analysis for biological datasets
* Automated testing with `testthat`
* Dynamic report generation using R Markdown

Each assignment mirrors real-world bioinformatics workflows where analysis logic, validation, and reporting are cleanly separated.

🔗 **Course Overview:**
[https://www.bu.edu/academics/cds/courses/cds-bf-591/](https://www.bu.edu/academics/cds/courses/cds-bf-591/)

---

## Repository Structure

Each assignment follows the same modular structure:

```
├── reference_report.html
├── main.R
├── README.md
├── report.Rmd
└── test_main.R
```

This structure enforces separation of concerns and reproducibility.

---

## Assignment Workflow

Each assignment follows a consistent analytical pipeline:

1. **Implement functions in `main.R`**

   * Write modular, reusable functions
   * Follow documented specifications

2. **Validate with `test_main.R`**

   * Source `main.R`
   * Run automated tests
   * Ensure correctness of outputs

3. **Generate report via `report.Rmd`**

   * Source implemented functions
   * Produce figures, summaries, and analyses
   * Knit to HTML output

4. **Compare to `reference_report.html`**

   * Ensure analytical consistency
   * Validate structure and expected outputs

This mirrors best practices used in research labs and production bioinformatics pipelines.

---

## Core Skills Developed

Through these assignments, students develop proficiency in:

* Writing clean, modular R functions
* Implementing statistical workflows
* Data visualization for biological interpretation
* Differential expression analysis foundations
* Tidy data principles
* Functional abstraction
* Debugging and test-driven development
* Reproducible report generation

---

## Reproducibility & Testing

Each assignment enforces:

* ✔ Encapsulation of logic in standalone functions
* ✔ Automated validation using test scripts
* ✔ Separation of computation from presentation
* ✔ Fully reproducible HTML reports
* ✔ No manual editing of analytical outputs

This structure reflects real-world computational biology workflows, where reproducibility and validation are critical.

---

## Installation & Setup

Clone the repository:

```bash
git clone https://github.com/your-username/bf591-r-biological-sciences.git
```

Install required R packages:

```r
install.packages(c(
  "tidyverse",
  "rmarkdown",
  "testthat"
))
```

To run an assignment:

```r
source("main.R")
source("test_main.R")
```

Then open and knit:

```
report.Rmd
```

---

## Academic Context

Course: **BF591 – R for Biological Sciences**
Institution: Boston University
Program: MS in Bioinformatics

This coursework provided foundational training in structured R programming and reproducible analytical workflows, forming the basis for advanced multi-omics and statistical genomics projects.

---

## 👩‍💻 Maintainer

**Neha Rao**
MS in Bioinformatics – Boston University

Focus Areas:

* Reproducible omics analysis
* Statistical genomics
* Transcriptomics workflows
* Computational biology pipelines
