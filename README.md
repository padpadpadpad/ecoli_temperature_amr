# Antimicrobial resistance in a warmer world: selection and transfer of kanamycin resistance is temperature dependent in *E. coli*

This repository contains raw phenotypic data, processing and analysis scripts, for the study on the selection and spread of kanamycin resistance in *E. coli* MG1655 and plasmids R1 and RP4.

- Initial submission archived on **Zenodo**: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23389520.svg)](https://doi.org/10.5281/zenodo.23389520)
- Preprint on [**bioRxiv**](https://doi.org/10.64898/2026.10.16.759609).

## Contents

- [Licensing](#licensing)  
- [Using this repository](#using-this-repository)  
- [Downloading the datasets](#downloading-the-datasets)  
- [Repository structure](#repository-structure)  
- [Scripts](#scripts)  
- [Data dictionary](#data-dictionary)

## Authors

- Zoltan Erdos
- Daniel Padfield  
- Tatiana Dimitriu  
- Ben Raymond  
- Angus Buckling

## Licensing

The code used for the analysis is licensed under an MIT license and has been developed with the help of GitHub Copilot, and potentially other LLMs. However, all output has been checked and modified by a human and we take full responsibility for any errors.

## Using this repository

This repository should be easy to run on any machine using RStudio projects or Positron folders. Both automatically set the root directory to the directory in which the project resides, so all analyses should run without altering paths.

- **RStudio:** create an `.Rproj` file in the base folder of this repository and open it as a project.  
- **Positron:** use *File \-\> Open Folder* to set the root directory to the base folder of the repository.

All packages used by a script are loaded at the top of each script, alongside a description of what each script does. RStudio and Positron usually offer to install missing packages; if not, the first code block of each script installs all necessary packages. All processed datasets are saved in `data/processed`, so each script can be run on its own.

## Downloading the saved models

GitHub is not designed to store large data files, so all saved model fits are stored in a snapshot of this repository on Zenodo. To download them and place them in `data/models` see **00\_download\_models.R** (Only available on the GitHub repository). All other data are included in this repository.

## Repository structure

| Folder | Contents |
| :---- | :---- |
| `data/` | All raw and processed datasets. See the [Data dictionary](#data-dictionary) for details of individual files. |
| `scripts/` | All code needed to recreate the processing and analyses. See [Scripts](#scripts) for details of each script. |
| `plots/` | All figures created for the manuscript, plus some extra figures that are not included. Main manuscript figures are named `figure_X.png` etc.; other figures have descriptive names. |
| `tables/` | Tables created for the manuscript. |

&nbsp;

## Scripts

Each script and a high-level overview of its function are listed below. For scripts that fit models, we provide code to load the fitted model saved in `data/models`, although you can refit the model if you wish.

> **Note:** Scripts numbered 01-05 must be run in order as they rely on sequential outputs. Run them separately on the dose-response dataset (`data/raw/growth_rates/DR*`) and the transconjugant dataset (`data/raw/growth_rates/T*`). The transconjugant growth data were obtained without antibiotics and no wells dried out, so **05\_clean\_data.R** does not produce any output and will not run for that dataset.

| Script | Description |
| :---- | :---- |
| **00\_flextpc\_check.R** | Compares the modified flexTPC model used throughout our analyses with the original formulation by Mauricio Cruz Loya ([flexTPC](https://mcruzloya.github.io/flexTPC/)). |
| **00\_functions.R** | Functions used in other analyses. |
| **01\_read\_od.R** | Lists and reads in plate reader output spreadsheets. |
| **02\_visualise\_runs.R** | Plots raw OD readings for each well in each plate. |
| **03\_process\_od.R** | Corrects the collated data for blanks, removes blanks and filters out biologically implausible data. |
| **04\_analyse\_gcplyr.R** | Removes noise at very low OD, then calculates and visualises growth metrics. |
| **05\_clean\_data.R** | Removes dried-out wells, *de novo* mutants etc. See the notes in the script and the methods section of the manuscript for details. |
| **brms\_tpc.R** | Bayesian analysis of the thermal performance curves of strains in the absence of antibiotics. |
| **brms\_plasmid.R** | Bayesian analysis of the thermal performance curves of plasmid transfer. |
| **brms\_drc.R** | Bayesian mixed-effects dose-response analysis of kanamycin growth data. Fits the model, calculates traits (EC50, EC95/MIC, NSEC, sub-MIC selection window) and selection estimates, and creates Figures 2 and 3\. |
| **invasion.R** | Bayesian analysis of the invasion assays. Creates Figure 4\. |
| **create\_figure\_1.R** | Recreates Figure 1 by combining outputs from **brms\_tpc.R** and **brms\_plasmid.R**. Creates table_tpc.png |

&nbsp;

## Data dictionary

This section describes the datasets in this repository.

### Raw growth rate measurements: `data/raw/growth_rates`

The files in this folder are output files from Multiskan SkyHigh microplate spectrophotometers (Thermo Fisher Scientific). Each file is a kinetic read of one 96-well plate. File names follow the pattern `DR` (dose response), then the date (day and month), temperature and replicate number, separated by underscores. For example, `DR1306_37_1` is a dose-response run on 13 June at 37 °C, replicate 1\.

Each `.xlsx` file has five sheets:

| Sheet | Description |
| :---- | :---- |
| Absorbance results | Absorbance of each well at each reading timepoint. |
| Raw data | The same data as Absorbance results, in long format, with measurement time instead of reading number. |
| General information | Run settings, including the serial number and firmware version of the machine. |
| Plate layout | Layout of the 96-well plate showing the strain in each well. Each column contains a different antibiotic concentration, from 64 μg/ml in column 1 down to 0 in column 11\. Column 12 is a blank. |
| Run log | Temperature record throughout the kinetic run. |

&nbsp;

The plate layouts have abbreviations to explain each well. The dose-response plate layout abbreviations (files beginning with **DR**) are as follows:

| Abbreviation | Strain |
| :---- | :---- |
| wt | Wild type MG1655 |
| R1 | Wild type \+ R1 |
| RP4 | Wild type \+ RP4 |
| nal.5 | nal.5 mutant |
| wt.c | Wild type without antibiotic |
| R1.c | Wild type \+ R1 without antibiotic |
| RP4.c | Wild type \+ RP4 without antibiotic |
| nal.5.c | nal.5 mutant without antibiotic |
| B | Sterile LB broth |

&nbsp;

Spreadsheets starting with **T** contain growth rates of transconjugants of R1 and RP4. These are used only to account for differences in strain-specific growth rate when estimating the plasmid transfer rate with the approximate extended Simonsen model. Some plates have additional controls, which can be seen in the Plate layout sheet. The transconjugant plate layout abbreviations are as follows:

| Abbreviation | Strain |
| :---- | :---- |
| wt | Wild type MG1655 |
| wt.r1 | Wild type \+ R1 |
| wt.rp4 | Wild type \+ RP4 |
| nal.5 | nal.5 mutant |
| R1 | R1 transconjugant |
| RP4 | RP4 transconjugant |
| B | Sterile LB broth |
| X | Uninoculated row |

&nbsp;

### Processed growth rate data: `05_kanamycin_metrics_analysis.csv`

Files in `data/processed` are intermediate outputs of the processing pipeline described above. The main analysis file, `05_kanamycin_metrics_analysis.csv`, is documented below.

| Column | Description |
| :---- | :---- |
| `file` | Name of the data file the record is from |
| `well` | Well the metrics correspond to |
| `serial_no` | Serial number of the OD reader |
| `id` | ID of the measured strain, taken directly from the OD plate layout (see above) |
| `temp` | (Constant) temperature of the kinetic OD read |
| `run` | Temperature (plate) replicate number |
| `threshold_density` | To remove noise at low ODs we set an arbitrary, very low threshold of 0.009 OD. We calculate the proportion of points in each well above the threshold and filter on this value. Low ODs are noisy but sometimes real (e.g. no growth ever), so wells are treated differently. `mostly_above`: at least 50% of points are above the threshold, so points below it were removed. `mostly_below`: fewer than 50% of points are above the threshold, so all data were kept |
| `Min_dens` | Minimum density |
| `lag_time` | Lag time |
| `gr` | Growth rate |
| `gr_time` | Per-capita growth rate over time |
| `gr_dens` | Per-capita growth rate over density |
| `max_dens` | Maximum density |
| `max_time` | Time to reach maximum density |
| `doub_time` | Doubling time |
| `final_time` | Final measurement time point |
| `final_dens` | Final density |
| `auc` | Area under the curve |
| `auc_48` | Area under the curve limited to the first 48 hours |
| `n` | Unique number identifier per measurement |
| `classification` | Programmatic classification of wells based on set criteria. Wells that dried out or showed signs of *de novo* evolution of resistance over the kinetic read are filtered out. The criteria are detailed in **05\_clean\_data.R** and in the methods section of the manuscript |
| `id2` | ID of the strain measured (`id` column), shortened from the plate layout notation. For example, for `wt.c` the `.c` control notation was removed because the antibiotic concentration is now attached (`conc_µl`) |
| `classification2` | Manual classification, based on our knowledge of the system, of wells we judged to have dried out or shown *de novo* evolution |
| `zoli_remove` | Wells to be removed, based on a comparison of `classification` and `classification2`. The final dataset should only contain the value `fine` |
| `well_num` | The number from the well name |
| `conc_µl` | Concentration of antibiotic in the well |

&nbsp;

Other key files are `wells_to_discard.xlsx` and `05_wells_to_check.csv`, explained below.

`wells_to_discard.xlsx`
This file contains wells that were identified as either dried out or show signs of *de novo* evolution. 

| Column | Description |
| :---- | :---- |
| `Plate_id` | Id of plate including day/month |
| `temp` | temperature of growth dose-response assay |
| `run` | Length of the plasmid transfer assay (hours) |
| `rep` | Temperature of the plasmid transfer assay (°C) |
| `comment` | Any notes prior to analysis |

`05_wells_to_check.csv`
Wells to check DP (Daniel Padfield) and ZE (Zoltan Erdos) contains the few wells of the dose repsonse assays where our manual checks didnt agree with our programmatic checks. We use our expertise of the system to classify these independently than compare the outcome. Wells were only removed if we independently agreed that they were failed (due to drying out or de novo evolution).

| Column | Description |
| :---- | :---- |
| `File` | Id of plate including day/month |
| `well` | well under scrutiny |
| `classification` | Manual classification, based on our knowledge of the system, of wells we judged to have dried out or shown *de novo* evolution |

### Plasmid transfer: `conjugation.xlsx`

Count data from plasmid transfer assays obtained by dilution plating.

| Column | Description |
| :---- | :---- |
| `strain` | Plasmid used in the transfer assays: R1 or RP4 |
| `plate` | `DT`: selective plate for donors and transconjugants (50 μg/mL kanamycin). `R`: selective plate for recipients (64 μg/mL nalidixic acid). `T`: double selective plate for transconjugants (50 μg/mL kanamycin and 64 μg/mL nalidixic acid) |
| `dilution` | Serial 10-fold dilution at which countable colonies were found |
| `rep` | Biological replicate number |
| `time_h` | Length of the plasmid transfer assay (hours) |
| `temp` | Temperature of the plasmid transfer assay (°C) |
| `count` | Colony count (6 technical replicates per biological replicate) |
| `density` | Density in CFU/mL, calculated as CFU / (plated volume (10 μL) × dilution factor) |
| `date` | Date of the assay |

&nbsp;

### Selection for resistance (invasion from rare): `invasion.xlsx`

Count data from the invasion experiment obtained by dilution plating.

| Column | Description |
| :---- | :---- |
| `strain` | Plasmid used in the invasion experiment |
| `conc` | Antibiotic concentration in the wells |
| `plate` | `DT`: selective plate for donors and transconjugants (50 μg/ml kanamycin). `R`: selective plate for recipients (64 μg/ml nalidixic acid). `T`: double selective plate for transconjugants (50 μg/ml kanamycin and 64 μg/ml nalidixic acid) |
| `temp` | Temperature of the invasion experiment (°C) |
| `dil` | Serial 10-fold dilution at which countable colonies were found |
| `count` | Colony count (6 technical replicates per biological replicate at each temperature × concentration combination) |
| `density` | Density in CFU/ml, calculated as CFU / (plated volume (10 μl) × dilution factor) |
| `time_h` | Length of the invasion assay (hours) |
| `rep` | Biological replicate number |

&nbsp;

### Model fits: `data/models`

This folder contains the Bayesian models fitted for this project. They can be read in to create credible intervals of key parameters and figures without rerunning the fitting process.

| File | Description |
| :---- | :---- |
| **brms\_flextpc.rds** | Thermal performance curves for all strains in the absence of antibiotics. |
| **brms\_flextpc\_r1.rds** | Thermal performance curve for the R1 plasmid transfer assays. |
| **brms\_flextpc\_rp4.rds** | Thermal performance curve for the RP4 plasmid transfer assays. This model has not converged, likely because the data are poor. |
| **brms\_drc.rds** | Dose-response curves for wild-type *E. coli*, nalR, wt \+ R1 and wt \+ RP4. |
| **brms\_kana\_traits.rds** | Fit to 500 draws of the estimates of MIC, EC50, NSEC and the sub-MIC selection window calculated from **brms\_drc.rds**. As it is essentially 500 brms fits, it takes a long time to run and is a large file. |
| **brms\_kana\_sc.rds** | Fit to 500 draws of the selection coefficients calculated from **brms\_drc.rds**. As it is essentially 500 brms fits, it takes a long time to run and is a large file. |
| **brms\_kana\_src.rds** | Fit to 500 draws of the selection rate coefficients calculated from **brms\_drc.rds**. As it is essentially 500 brms fits, it takes a long time to run and is a large file. |
| **beta\_invasion.rds** | How the proportion of resistant bacteria changes across temperatures and kanamycin concentrations in invasion assays. |
| **negbinom\_invasion.rds** | How the number of transconjugants changes across temperatures and kanamycin concentrations in invasion assays. |
| **negbinom\_invasion2.rds** | How the number of transconjugants changes across temperatures and kanamycin concentrations in invasion assays. Testing a simpler model structure. |
| **negbinom\_invasion3.rds** | How the number of transconjugants changes across temperatures and kanamycin concentrations in invasion assays. Testing a simpler model structure. |

&nbsp;