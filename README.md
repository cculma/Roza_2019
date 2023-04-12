# Roza2019

This project summarize the data analysis for Roza2019 field.

## Phenotypic data

Roza2019 is a half-sib population with 436 alfalfa half-sib families each one contains at least 18 individuals which were developed by nested associated mapping (NAM).

Roza 2019 NAM population was constructed using a drought resistant parent population from WA and UT (1608) and a drought sensitive parent population from Guardsman II synthetic cultivar (1613). Guardsman II was planted in the greenhouse and after an evaluation protocol and selected 30 genotypes that were most susceptible to drought stress. Seed of F1 population (1622) were collected of 1613 parent and labeled A-F while seeds collected of 1608 parent were labeled G-L.

## Yield

Harvest has been done during 2020, 2021, 2022 and posible 2023. Data modeling can be done using two stages approach or random regression model. Because this is an experiment with repeated measurements or longitudinal data, Random Regression Analysis is a good option.

|          | Single | Trial |       | Stagewise |          |             |
|----------|--------|-------|-------|-----------|----------|-------------|
| Roza2019 | 2020   | 2021  | 2022  | gen_month | gen_year | gen_overall |
|          | May    | May   | -     | ST2_may   | ST3_20   | ST4_Yi      |
|          | Jun    | Jun   | Jun   | ST2_jun   | ST3_21   | -           |
|          | Jul    | Jul   | Jul   | ST2_jul   | ST3_22   | -           |
|          | Aug    | Aug   | Aug   | ST2_aug   | -        | -           |
|          | Sep    | Sep   | Sep   | ST2_sep   | -        | -           |
|          | Total  | Total | Total | -         | -        | -           |

## R scripts

1. Roza2019_GWAS.R is the main file wich runs GWASpoly analysis using a phenotypic matrix of $400 \times 23$ and a genotypic matrix of $82,156 \times 501$ in GWASpoly format 'ACGT'.

2. Install and run `updog` to check LD decay.


