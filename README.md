# Roza2019

This project summarize the data analysis for Roza2019 field.

## Phenotypic data

Roza2019 is a half-sib population with 436 alfalfa half-sib families each one contains at least 18 individuals which were developed by nested associated mapping (NAM).

Insert table of population.

Roza 2019 NAM population was constructed using a drought resistant parent population from WA and UT (1608) and a drought sensitive parent population from Guardsman II synthetic cultivar (1613). Guardsman II was planted in the greenhouse and after an evaluation protocol and selected 30 genotypes that were most susceptible to drought stress. Seed of F1 population (1622) were collected of 1613 parent and labeled A-F while seeds collected of 1608 parent were labeled G-L.

## Yield

Harvest has been done during 2020, 2021, 2022 and possibly 2023. Data modeling can be done using two stages approach or random regression model. Because this is an experiment with repeated measurements or longitudinal data, Random Regression Analysis is a good option.

|          | Single | Trial |       | Stagewise |          |             |
|----------|--------|-------|-------|-----------|----------|-------------|
| Roza2019 | 2020   | 2021  | 2022  | gen_month | gen_year | gen_overall |
|          | May    | May   | -     | ST2_may   | ST3_20   | ST4_Yi      |
|          | Jun    | Jun   | Jun   | ST2_jun   | ST3_21   | -           |
|          | Jul    | Jul   | Jul   | ST2_jul   | ST3_22   | -           |
|          | Aug    | Aug   | Aug   | ST2_aug   | -        | -           |
|          | Sep    | Sep   | Sep   | ST2_sep   | -        | -           |
|          | Total  | Total | Total | -         | -        | -           |

Note: obtain percentage of each yield in total yield.
Harvest 1 + Harvest 2 + Harvest 4 = 100% (Season total yield)

## R scripts

1. Roza2019_GWAS.R is the main file wich runs GWASpoly analysis using a phenotypic matrix of $400 \times 23$ and a genotypic matrix of $82,156 \times 501$ in GWASpoly format 'ACGT'.

2. Install and run `updog` to check LD decay.

3. Try SNP-based association to Gene-based association.

## Definitions

1. SNV: Single nucleotide variant is a variation of a single nucleotide in a population's genome.
2. SNP: Single nucleotide polimorphism is the same variation of a single nucleotide in a population's genome, but it is limeted to germline DNA and must be present in at leat 1% of the population. MAF filtering 0.01.

## Generation of VCF

Please check file `bash_scripts/4.1_ngsep_NMGS.sh`

1. Generate the raw VCF file with the function `MultisampleVariantsDetector`

- `java -Xmx50g -Xms40g -jar ${NGSEP} MultisampleVariantsDetector -maxAlnsPerStartPos 100 -maxBaseQS 30 -ploidy 4 -psp -knownSTRs ${STR} -r ${GENOME} -o Roza2019_01.vcf `

2. Filter the raw VCF file using multiple filtering parameters with the function `VCFFilter`

- `java -Xmx50g -Xms45g -jar ${NGSEP} VCFFilter -q 40 -s -fi -m 250 -minRD 8 -i Roza2019_01.vcf -o ../4_vcf_monoploid/Roza2019_02.vcf`

Oct 11, 2023 2:58:08 PM ngsep.vcf. \
VCFFilter logParameters \
INFO: Input file: Roza2019_01.vcf \
Output file: ../4_vcf_monoploid/ \ Roza2019_02.vcf \
Genotype filters \
Minimum genotype quality: 40 `-q 40` \
Minimum read depth: 8 `-minRD 8`\
Variant context filters \
Population data filters \
Minimum samples genotyped: 250 `-m 250` \
Keep only biallelic SNVs `-s` \
Filter sites where only one allele is observed in the population `-fi`

3. GWASPoly format

- `java -Xmx50g -Xms45g -jar ${NGSEP} VCFConverter -GWASPoly -i Roza2019_03_imputed.vcf -o Roza2019_04`
