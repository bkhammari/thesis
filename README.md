# How Equal is Equal Opportunity? 🎓🇧🇷
## Gender, Ethnic and Regional Disparities in Labor Market Outcomes of Higher Education Policies in Brazil

Author: Baha Eddine Khammari
M.Sc. Economics, University of Cologne

This repository contains the code, documentation, and thesis materials for my Master's Thesis in Economics at the University of Cologne. The project studies how Brazil’s **ProUni** (*Programa Universidade para Todos*) affected early-career labour-market outcomes and whether these effects differed by gender, demographic background, and geography.

## Overview

The thesis examines the labour-market consequences of higher-education expansion in Brazil using a microregion-by-year panel covering **558 Brazilian microregions** over **2005–2019**. The empirical design combines ProUni administrative records, RAIS formal labour-market microdata, and 2010 Census baseline information to estimate how variation in local scholarship intensity translated into differences in wages, employment, and early-career labour-market integration.

Substantively, the project asks whether expanding access to higher education produced equal gains across groups, or whether existing structural inequalities in Brazil’s labour market shaped who benefited most. The focus is therefore not only on average treatment effects, but also on distributional heterogeneity across gender, ethnicity, and regional labour markets.

## Research Questions

- What is the effect of ProUni expansion on formal-sector earnings and employment?
- Do labour-market effects differ across low- and high-exposure local labour markets?
- Are the gains from higher-education expansion distributed differently across:
  - Gender
  - Demographic background / ethnicity
  - Geographic regions
- Do the results point to scarcity rents, saturation effects, or broader general-equilibrium spillovers?

## Empirical Strategy

The main analysis uses modern staggered-adoption Difference-in-Differences methods based on **Callaway & Sant’Anna (2021)**. Treatment is constructed as a lagged cumulative measure of ProUni scholarship density at the microregion level, and the baseline empirical design groups microregions into dose strata to estimate interpretable event-study effects.

The repository also includes:
- Doubly robust specifications with pre-treatment covariates
- Cohort-based designs comparing younger and older workers
- TWFE benchmark models using `fixest`
- Continuous-dose robustness exercises
- Sensitivity checks for lag structure, cutoff definitions, and potential spillovers

## Data

### Primary Sources
- **ProUni administrative records**: scholarship allocations by location and year
- **RAIS** (*Relação Anual de Informações Sociais*): matched employer-employee administrative data on Brazil’s formal labour market
- **IBGE Census 2010**: demographic baseline variables used for normalization and covariates
- **Geographic concordances** via `geobr` and Brazilian territorial directories

### Access and Processing
Data extraction is handled in **R** through [`basedosdados`](https://basedosdados.org/), which allows direct querying from Google BigQuery. The workflow harmonizes municipality-level information into **microregion-level local labour markets**, which form the main unit of analysis in the thesis.

## Main Outcomes

The repository contains code to construct and analyse:
- Average log formal wages
- Formal employment counts / extensive-margin outcomes
- Age-cohort wage gaps
- College-share first-stage outcomes
- Heterogeneity estimates by gender and ethnicity
- Regional comparisons and robustness exercises

## Tech Stack

- **Language**: R
- **Econometrics**: `did`, `fixest`, `modelsummary`
- **Data work**: `tidyverse`, `data.table`, `basedosdados`
- **Geography**: `geobr`
- **Visualisation**: `ggplot2`
- **Writing**: LaTeX / Overleaf

## Reproducibility

The project is designed around a reproducible workflow:
1. Query and harmonise raw administrative data
2. Construct treatment intensity and outcome panels
3. Run baseline and heterogeneity models
4. Export tables, figures, and audit logs for the thesis appendix

The current pipeline writes outputs to folders such as:
- `FiguresFinal/`
- `TablesFinal/`
- `DocumentationFinal/`

Because some inputs depend on external data access credentials and BigQuery permissions, full replication may require a valid `basedosdados` billing project and access to the underlying public datasets. In this thesis, I used `microdadosbrazil` as a billing ID for the Google Cloud Queries.

