# Artificial Intelligence based Histological Markers Predict Response to Neoadjuvant Chemotherapy in Triple-Negative and HER2+ Breast Cancers
by
Rawan Albusayli,
Liselore M. Janssen,
etc

## Abstract
Neoadjuvant Chemotherapy (NAC) is a commonly employed treatment strategy for patients with Triple-Negative Breast Cancer (TNBC) and HER2+ Breast Cancer (BC) cases, aimed at reducing tumor size before surgery. However, not all patients respond favorably to NAC and some experience adverse outcomes. Therefore, identifying predictive factors that ascertain treatment benefits becomes crucial in steering clinical decisions, ideally before the commencement of the NAC treatment. In this study, we explore the correlation of novel deep learning-based histological markers extracted from digital images of routinely stained slides of pre-treatment biopsies of 100 patients with NAC response in TNBC and HER2+ subtypes. We find that a higher value of our artificial intelligence (AI)-derived quantitative measure of tumor associated stroma (Digi-TAS score) is associated with less favorable responses to NAC (RCB > 0), with Area Under the ROC Curve (AUROC) value of 0.72 for the overall cohort and 0.81 for the TNBC cohort. Conversely, a higher value of our score for stromal tumor-infiltrating lymphocytes (or the Digi-sTILs score) is associated with more positive outcomes, especially in cases of HER2+ BC, with an AUROC value of 0.65. Our innovative AI approach offers promise for identifying patients likely to respond to NAC and could potentially assist clinicians in treatment decisions in future.

## Software implementation

The source code used to generate the results in the paper is provided as separate files for each process, which can be run separately to create graphs, calculate scores, and finally make the statistical calculations with NAC response metrics.

## Getting the code

You can download a copy of all the files in this repository by cloning the
[git](https://git-scm.com/) repository:

    git clone https://github.com/Rawanrb/NAC_Paper.git

## Dependencies

You'll need a working Python environment to run the code.
The recommended way to set up your environment is through the
[Anaconda Python distribution](https://www.anaconda.com/download/) which
provides the `conda` package manager.
Anaconda can be installed in your user directory and does not interfere with
the system Python installation.
The required dependencies are specified in the file `environment.yml`.

We use `conda` virtual environments to manage the project dependencies in
isolation.
Thus, you can install our dependencies without causing conflicts with your
setup (even with different Python versions).

Run the following command in the repository folder (where `environment.yml`
is located) to create a separate environment and install all required
dependencies in it:

    conda env create


## Reproducing the results

Before running any code you must activate the conda environment:

    source activate ENVIRONMENT_NAME

or, if you're on Windows:

    activate ENVIRONMENT_NAME

This will enable the environment for your current terminal session.
Any subsequent commands will use software that is installed in the environment.
To execute and test the code, uncomment the lines in `code/runme.R` or `code/runme.py`.

## License 

All source code is made available under a BSD 3-clause license. You can freely
use and modify the code, without warranty, so long as you provide attribution
to the authors. See `LICENSE.md` for the full license text.

The manuscript text is not open source. The authors reserve the right to the
article content.
