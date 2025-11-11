## **README — biniLasso: Automated Cut-point Detection via Cumulative Binarization**

### **Overview**

This repository provides the R code, datasets, and vignette examples for the **biniLasso** and **miniLasso** methods, introduced by Safari, Helisaz, and Loewen (2025). These methods are designed for **automatic cut-point detection** in high-dimensional survival data, offering interpretable and efficient alternatives to existing approaches such as **Binacox**. Both methods are implemented within the Cox proportional hazards framework and make use of **cumulative binarization** combined with **L1-penalized regression** to detect optimal thresholds for continuous features. The sparse variant, **miniLasso**, integrates **uniLasso** regularization to enhance sparsity while preserving univariate effect directions.

---

### **Purpose**

The vignette demonstrates how to:

1. Load and prepare high-dimensional cancer datasets (GBM, BRCA, and KIRC).
2. Generate candidate cut-points using the `cumBinarizer()` function.
3. Identify optimal thresholds with `opt_cuts_finder()` under both biniLasso and miniLasso.
4. Compare the performance of fitted Cox models using AIC, IBS, and C-index metrics.
5. Restrict the number of thresholds per variable using `opt_fixed_nCuts()` for better interpretability.

Each step in the vignette is fully reproducible, and all functions are part of the **biniLasso R package**.

---

### **Data**

The vignette uses three preprocessed TCGA datasets:

* `gbm_fnl.rds` — Glioblastoma multiforme (GBM)
* `brca_fnl.rds` — Breast invasive carcinoma (BRCA)
* `kirc_fnl.rds` — Kidney renal clear cell carcinoma (KIRC)

These datasets contain gene expression measurements, survival times, and censoring information. The code assumes these files are available in the `data/` directory.

---

### **Code Structure**

* **`biniLasso vignette.Rmd`** — Main vignette file illustrating the entire workflow.
* **`data/`** — Folder containing processed TCGA datasets and Binacox comparison results.
* **`R/`** — Source code for biniLasso and related helper functions.
* **`results/`** — Output tables showing detected cut-points and performance metrics.

---

### **Notes on Fixed Number of Thresholds**

When the maximum number of cut-points per gene is fixed (for example, at two), there are cases where a variable ends up with only one threshold or an `NA`. This happens when the model detects overlapping or nearly identical cut-points, which are merged into a single effective boundary, or when no meaningful threshold is found after penalization. In the latter case, the variable’s effect is likely linear or weakly associated with the outcome, so no threshold is selected. These results are expected and reflect the model’s regularization process, which automatically removes redundant or non-informative cut-points to ensure simplicity and interpretability.

---

### **Reference**

Safari, A., Helisaz, H., & Loewen, P. (2025). *biniLasso: Automated cut-point detection via sparse cumulative binarization.* arXiv preprint [arXiv:2503.16687](https://arxiv.org/abs/2503.16687).

---

### **Contact**

For questions, suggestions, or bug reports, please contact:
**Abdollah Safari** – [a.safari@ut.ac.ir](mailto:a.safari@ut.ac.ir)
**Hamed Helisaz** – [hamed.helisaz@ubc.ca](mailto:hamed.helisaz@ubc.ca)

---

Would you like me to tailor this README for inclusion on GitHub (with Markdown formatting, links, and badges), or keep it as a plain text document for internal use (e.g., within the vignette folder)?
