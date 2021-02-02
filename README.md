# frictions-JHR

This repository contains code to replicate the analysis and results in "Labor Market Frictions and Moving Costs of the Employed and Unemployed" published in the *Journal of Human Resources*

## Data access requirements
Access to a Research Data Center (RDC) of the U.S. Census Bureau is required to fully replicate the analysis. I include all code needed to replicate the analysis. Note that file paths have been redacted in keeping with Census Bureau disclosure policies.

## Software requirements
The data and analysis makes use of the following software packages:

- SAS
- Stata
- R
- Matlab
- Julia

## Data Sources
The primary data source is a confidential version of the 2004 and 2008 panels of the Survey of Income and Program Participation (SIPP).

This data source is enriched with the following other data sources:

- Proprietary data from the American Chamber of Commerce Researchers Association (ACCRA) on cost of living in various U.S. Cities. 
    * I obtained this data from Christopher Timmins (Duke Univ.) and am not able to share it publicly due to a data agreement.
- Publicly available data from the American Community Survey (ACS) on rental rates across various locations in the U.S.
    * I access the ACS data from the Integrated Public Use Microdata Series (IPUMS). I am not able to share the ACS microdata publicly due to the terms of use with IPUMS.
- Publicly available data on county unemployment rates over time, obtained from the Bureau of Labor Statistics Local Area Unemployment Statistics (BLS-LAUS)
- These data sources are used to create a county-level data set which is locatd in `Data/Input/AccraData/finaldata/county_city_chars.dta`

Once the data for these input sources is created, I merge it with the confidential SIPP data.

For descriptive statistics in the paper, I also make use of publicly available SIPP data, obtained from the National Bureau of Economic Research (NBER)

## Code for replication
Because the code uses five different languages, there is no single "master" script that can be executed to completely replicate every result. Therefore, I walk through the steps for data creation and estimation:

### Custom code files
I store custom Linux commands and custom Matlab commands in the `Functions/` folder.

### Input data
Directions for creating the county-level data are located in `Data/Input/README.md`

### Public-use SIPP data
Data was downloaded from NBER, see `Data/PublicSIPP/README.md` for complete details. I exclude the binary data files in order to minimize disk space. These binary SIPP data files are also available [here](https://github.com/tyleransom/EER_CWP/tree/master/)

### Assembling restricted-use SIPP data
Run `Data/DataScriptsRDC/data_create.bash`. Filepaths are redacted due to Census Bureau disclosure requirements.

### Replicating Figures and Tables presented in the paper
#### Descriptive Statistics from Public-Use Data (Figure 1, Table 3; Figure A2, Table A6; Figure A3, Table A7; Table A8)
See code on lines 425-508 of `Data/PublicSIPP/Analysis/an_sipp_men.do` and similar lines in `Data/PublicSIPP/Analysis/an_sipp*.do`. Graphics are output as `.eps` files in the working directory, but are also copied over to `Graphics/` for simplicity.

#### Descriptive Statistics from Restricted-Use Data (Tables 1-2, A1-A2)
Sample selections are made in `Data/DataScriptsRDC/{cr_clean_SIPP.do,cr_clean_SIPP_2008.do}` and (manually) output to `Tables/SampleSelection.xls`.

#### Structural Estimates from Restricted-Use Data (Tables 4-6)
Run `Estimation/estimateFrictAR1IntBeta9_ub_only_1.m` which (manually) outputs to `Tables/structural_estimates.xls`.

#### Moving Cost Calculations (Table 7)
Run `Estimation/solveMC.jl`. Results print to console.

#### Model Fit (Tables 8-9)
Run `Estimation/modelFitFrictBeta9_ub_only_1.m` which (manually) outputs to `Tables/modelFit.xls`.

#### Counterfactuals (Figure 2 and Figures A6-A11)
To generate the data underlying the figures, run `Estimation/cflFrictAAmlmYYYYBornplttBeta9_ub_only_1.m` which (manually) outputs to `Tables/{allCflsLong.xls,Cfl_pool.xls,Counterfactuals.xls}`. To generate the figures, run `Estimation/Postestimation/GrapherCfls.do`.

#### Wage Elasticity to the Firm (Table A13)
Run `Estimation/concentration_markov_35_700.jl` and select rows of output from `Estimation/results_all_35_700_static.csv`.

#### Map (Figure A1)
Run `Graphics/Map/USA.R` which outputs `Graphics/map.eps`.

#### Location Characteristics (Tables A9-A12)
Run `Estimation/backOutCorrelationsFrictBeta9.m` and `Estimation/directionRegressionConcatenator.m` which (manually) outputs in `Tables/CityCorrelations.xls`.


