# Antimicrobial resistance patterns and multidrug resistance in a laboratory isolate dataset

> **Data warning.** The dataset used here is `example_isolates` from the R package
> `AMR`. The package describes it as *randomised, fictitious data that reflects reality and can
> be used to practise AMR data analysis*. It is **not** a real clinical sample collection.
> The pipeline, statistics and code are real and reproducible. The biological conclusions are
> **not** valid for any real population. Some results also break known microbiology (see
> Limitations). To apply this project to real data, see *Using real data* below.

## Abstract

Antimicrobial resistance (AMR) limits treatment options in hospitals, and laboratories track
resistance patterns and multidrug resistance (MDR) to guide empirical therapy. We analysed a
dataset of 2,000 isolates from 981 patients, tested against 40 antimicrobials (2002 to 2017). After
reducing to the first isolate per patient and species (1,281 isolates), we estimated resistance
prevalence for each drug, compared resistance by species, and defined MDR as resistance to three or
more antimicrobial classes. MDR was common (89% of isolates with enough data to classify). Adjusted
logistic regression showed no clear association with ward type or age. *S. aureus* isolates were
less likely to be MDR. A logistic regression predicted ciprofloxacin resistance from other antibiotic
results with 93% accuracy (AUC 0.90), against an 86% majority-class baseline, but it missed 37% of
resistant isolates. These are results for the dataset as provided, and they cannot be generalised to
a real hospital.

## 1. Introduction

Resistance to antimicrobials is a growing public-health problem. Clinicians often must start
treatment before susceptibility results return, so the pattern of resistance in a hospital guides the
choice of empirical antibiotics. Two standard questions are how common resistance is for each drug,
and how often isolates are multidrug resistant. MDR is usually defined as resistance to at least
three antimicrobial classes (Magiorakos et al., 2012).

Our aim was a data-driven question: **How common is resistance to each antimicrobial, is multidrug
resistance associated with ward type or patient characteristics, and can resistance to ciprofloxacin
be predicted from other antibiotic results?** This is a stewardship question, and the answer would
matter to a laboratory that wants to predict resistance without waiting for a full panel.

## 2. Data

**Source.** The `example_isolates` dataset in the `AMR` R package (version 3.x, GPL-2; Berends et al.,
2022). The file is copied into `data/` with the package's licence (`data/AMR_LICENSE.txt`). The package
also provides `antimicrobials`, used here to map each drug to its antimicrobial class.

**Size and variables.** 2,000 isolates from 981 patients, received between 2 January 2002 and
28 December 2017. Variables include the isolate `date`, `patient` ID, `age`, `gender`, `ward` (Clinical,
ICU or Outpatient), the microorganism code `mo`, and 40 antimicrobial columns coded S (susceptible),
I (intermediate) or R (resistant). Within the 40 columns, the number of isolates tested varies from 34
to 1,894 in the full dataset. There are no missing ages or genders in this copy.

**Quality.** Every antimicrobial value is S, I, R or missing. There are 552 repeated patient-species-date
combinations, so the first-isolate step below is important. The original package does not state
whether the randomisation preserved any real structure. We therefore cannot treat the statistical
associations as biology.

**Ethics and licence.** The data are fictitious, so no patient is identifiable. Use follows the GPL-2
licence of the `AMR` package.

## 3. Methods

The workflow is in `code/01_analysis.R`.

**Import and QC.** We loaded both datasets, checked that every antimicrobial code appears in the
reference table, checked the S/I/R coding, and counted missing values and duplicates.

**Cleaning: first isolate.** Repeat isolates from one patient with the same species are not
independent. We kept the earliest isolate for each patient and species, which is the standard
first-isolate approach. This gave 1,281 isolates. Species codes were converted to names for the six
most common species; the rest are labelled "Other species".

**Prevalence.** For each drug, the percentage of tested isolates classed R was computed with an exact
binomial 95% confidence interval. Drugs with 100 or more tests were plotted (Fig. 1). Resistance was
compared across species for the four most common species (Fig. 2).

**Multidrug resistance.** Each drug was assigned to its first listed antimicrobial class. An isolate was
classed MDR if it was resistant to at least three classes, provided at least three classes had been
tested. Isolates with fewer than three tested classes were left unclassified, so "not MDR" means
something.

**Association.** A logistic regression of MDR on ward, age, gender, and an *S. aureus* indicator
estimated adjusted odds ratios. Species was reduced to an *S. aureus* indicator because some species had
almost no non-MDR isolates, which caused complete separation in an earlier attempt.

**Prediction of ciprofloxacin (CIP) resistance.** The target was CIP resistance (R versus S or I).
Features were the R status of non-fluoroquinolone drugs tested on at least 1,000 isolates. Fluoroquinolones
were excluded so the model could not simply copy the target's own class. Only isolates with complete
feature data were used (396 isolates). The data were split 70/30 with stratification. A logistic
regression and an RBF-kernel support vector machine (SVM, tuned by 10-fold cross-validation on the training
set) were compared to a majority-class baseline.

**R packages.** `tidyverse` for data handling and graphics, `pheatmap` for the heatmap, `caret` for the split
and confusion matrix, `e1071` for the SVM. Versions are in `results/session_info.txt`.

## 4. Results

### 4.1 Data and first isolates

Of the 2,000 isolates, 1,281 remained after keeping the first isolate per patient and species. Across all
tested antimicrobials, 40.8% of results were R.

### 4.2 Resistance prevalence

Resistance was highest for kanamycin (100% of 328 tests), colistin (85.1%, n = 1,070), teicoplanin (72.5%)
and ceftazidime (71.9%). The full table with confidence intervals is in `results/drug_resistance_prevalence.csv`.

<img width="1050" height="1200" alt="fig1_prevalence" src="https://github.com/user-attachments/assets/25a3d299-dc70-4af6-8b71-3846d387b612" />


*Fig 1. Percentage resistant for each drug with at least 100 tests. Bars are 95% exact confidence intervals.*

Resistance differed by species (Fig. 2). Several species show resistance rates of 100% to
multiple drugs, which is unusually extreme; the Limitations section discusses this.

<img width="3000" height="1350" alt="fig2_species_heatmap" src="https://github.com/user-attachments/assets/1a764e00-14c8-422b-816e-92047f319b2f" />


*Fig 2. Percentage resistant for the four most common species against the 15 most-tested drugs. Cells with
fewer than 20 tests are grey.*

### 4.3 Multidrug resistance

Among 1,268 isolates with at least three antimicrobial classes tested, 89.2% were MDR. MDR was similar across
wards: Clinical 88.5% (687 of 776), ICU 90.0% (367 of 408) and Outpatient 91.7% (77 of 84) (Fig. 3).

<img width="900" height="750" alt="fig3_mdr_by_ward" src="https://github.com/user-attachments/assets/dc564ee0-6ef0-4b05-bddf-e4ee6ca7777a" />


*Fig 3. MDR prevalence by ward type, with approximate 95% confidence intervals.*

In the adjusted model (Table 1), neither ward type nor age was associated with MDR. *S. aureus* was
strongly less likely to be MDR than other species (odds ratio 0.14, 95% CI 0.09 to 0.21).

| Term | Adjusted odds ratio | 95% CI | p |
|---|---|---|---|
| Ward: Clinical (vs Outpatient) | 0.60 | 0.26 to 1.40 | 0.24 |
| Ward: ICU (vs Outpatient) | 0.64 | 0.27 to 1.55 | 0.33 |
| Age (per year) | 1.00 | 0.99 to 1.01 | 0.90 |
| Male (vs female) | 0.82 | 0.57 to 1.19 | 0.30 |
| *S. aureus* (vs other species) | 0.14 | 0.09 to 0.21 | < 0.001 |

*Table 1. Adjusted odds ratios for MDR (`results/mdr_logistic_odds_ratios.csv`).*

### 4.4 Predicting ciprofloxacin resistance

The model used 396 isolates with complete data, of which 54 (13.6%) were CIP-R. In the 118-isolate test set,
a majority-class baseline achieves 86.4% accuracy.

| Model | Accuracy | Sensitivity | Specificity | AUC |
|---|---|---|---|---|
| Logistic regression | 0.932 | 0.625 | 0.980 | 0.904 |
| RBF SVM (tuned) | 0.915 | 0.438 | 0.990 | 0.951 |

*Table 2. Test-set performance (`results/cip_model_performance.csv`). The test set has only about 16 resistant
isolates, so sensitivity estimates are imprecise.*

Both models beat the baseline on accuracy and AUC. Their sensitivity is low, so they miss many resistant isolates.
That matters clinically, because a missed resistant isolate could lead to an ineffective empirical choice.

## 5. Discussion

**Interpretation.** On this dataset, resistance is common for several drugs, and MDR is very common
overall. Ward type does not explain MDR. A ciprofloxacin model using other antibiotic results gives good
overall accuracy and AUC, but poor sensitivity for resistant isolates. The high MDR rate is a direct
consequence of the definition: with 40 drugs and 3 or more classes, most isolates qualify.

**Relation to the real-world problem.** A stewardship laboratory would use this type of analysis to
decide which drugs to report and which to test by default. The main finding for a real hospital would
depend on real data. This dataset cannot provide that answer.

**Limitations.**
1. *The data are fictitious.* The `AMR` package describes this dataset as randomised. Any real association
   in it may be chance, and the statistics answer only questions about this file.
2. *Some results are unusually extreme.* Several species show 100% resistance to multiple drugs, and
   kanamycin resistance is 100% across all 328 tests. Rates this uniform are rare in real surveillance data,
   so they should be checked against laboratory breakpoints before any interpretation. Some patterns are
   plausible; for example, the colistin resistance of *S. aureus* reflects intrinsic resistance of Gram-positive
   bacteria to colistin. These features point to synthetic data rather than real susceptibility.
3. *MDR definition.* MDR depends on the chosen class grouping and on how many drugs were tested. Using the first
   listed class of each drug is a simplification.
4. *Small test set for prediction.* The CIP test set has about 16 resistant isolates, so sensitivity has a wide
   margin of error. Results from one random split are not stable enough to rank the two models.
5. *Biological plausibility of features.* Features were chosen only by test coverage. A real model should use
   known resistance mechanisms and cross-resistance patterns, such as those in EUCAST expert rules.
6. *Missing and untested drugs.* Drugs are tested on different numbers of isolates, which means the
   complete-case analysis keeps only a subset of the data.

## 6. Conclusion and future work

In the dataset as supplied, resistance is widespread and MDR is common. Ward type and age are not associated
with MDR, and *S. aureus* has lower MDR odds. Ciprofloxacin resistance can be predicted with good overall
accuracy, but the model misses many resistant isolates.

Future work should:
- run the same pipeline on real laboratory surveillance data, such as WHONET exports or the EARS-Net
  surveillance data;
- use EUCAST or CLSI breakpoints and expert rules to check the S/I/R coding;
- add genotype features, such as resistance genes from AMRFinderPlus or CARD, for a genotype-phenotype model;
- use cross-validation repeated over many splits and report confidence intervals for sensitivity;
- weight the classifier towards sensitivity, because missing resistant isolates is the costly error.

## 7. Using real data

To apply this project to a real dataset, export isolates to a CSV with columns `patient`, `date`, `mo` (species),
`ward`, `age`, `gender` and one column per antimicrobial coded S, I or R. Then point `code/01_analysis.R` at that
file in place of the `load()` calls, and check the drug codes against a reference table. Real surveillance
data are restricted, so check the licence before publishing results from it.

## 8. Reproducibility

**Requirements.** R 4.3.3 or later. CRAN packages `tidyverse`, `pheatmap`, `caret`, `e1071`. Data are included
in `data/`, so no download is needed.

**Steps.** From the project root:

```bash
Rscript code/01_analysis.R
```

The script sets `set.seed(2026)` and uses relative paths. It writes tables to `results/`, figures to
`results/figures/`, and `results/session_info.txt` with package versions.

**Rendered output.** Plain R script; no HTML or PDF render is included.

## 9. References

Berends, M. S., Luz, C. F., Friedrich, A. W., Sinha, B. N. M., Albers, C. J., & Glasner, C. (2022). AMR: An R
package for working with antimicrobial resistance data. *Journal of Statistical Software*, 104(3), 1-31.
https://doi.org/10.18637/jss.v104.i03

Berends, M. S., et al. *AMR: Antimicrobial Resistance Data Analysis* (R package). The `example_isolates`
and `antimicrobials` datasets are distributed with this package under GPL-2.
https://github.com/msberends/AMR

Chang, C.-C., & Lin, C.-J. (2011). LIBSVM: A library for support vector machines. *ACM Transactions on
Intelligent Systems and Technology*, 2(3), 27.

Kolde, R. (2019). *pheatmap: Pretty Heatmaps* (R package version 1.0.12).
https://CRAN.R-project.org/package=pheatmap

Kuhn, M. (2008). Building predictive models in R using the caret package. *Journal of Statistical Software*,
28(5), 1-26.

Magiorakos, A.-P., et al. (2012). Multidrug-resistant, extensively drug-resistant and pandrug-resistant
bacteria: an international expert proposal for interim standard definitions for acquired resistance.
*Clinical Microbiology and Infection*, 18(3), 268-281.

Meyer, D., Dimitriadou, E., Hornik, K., Weingessel, A., & Leisch, F. (2023). *e1071: Misc functions of the
Department of Statistics, Probability Theory Group* (R package version 1.7-14).

R Core Team (2024). *R: A language and environment for statistical computing*. R Foundation for Statistical
Computing, Vienna, Austria. https://www.R-project.org/

Wickham, H., et al. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686.
