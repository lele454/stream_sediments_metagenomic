# Association between land-use intensity and redox reorganization in prokaryotic guilds of Atlantic Rainforest stream sediments

Code and derived tables for *Association between land-use intensity and redox reorganization in prokaryotic guilds of Atlantic Rainforest stream sediments
*.

Fourteen stream sediment metagenomes from the coastal plain of São Paulo,
Brazil: seven inside the Juréia-Itatins Ecological Station and seven in
catchments under human use. The analysis compares the two groups and relates
community composition to riparian land cover.

Everything here runs from the Bracken genus reports and the KEGG orthologue
table onward. Read processing, assembly, gene prediction and taxonomic
classification are described in the manuscript.

## Layout

```
scripts/   numbered pipeline, run in order
tables/    derived tables, as produced by the scripts
docs/      decisions that are not obvious from the code
```

## Requirements

```bash
conda create -n streams -c conda-forge \
    python=3.12 pandas numpy scipy scikit-learn matplotlib openpyxl
conda activate streams
```

## Inputs

Two files the scripts expect and this repository does not carry:

**Bracken genus reports**, one per sample, in a single directory. The sample
name must appear in the filename. Each is the tab-separated output of Bracken
at genus level, with the columns `name`, `taxonomy_lvl` and `new_est_reads`.

**KEGG orthologue table**, `.xlsx` or `.csv`, with KO identifiers in the first
column and one column per sample, in abundance per 100,000 annotated ORFs.

**Landscape table**, semicolon-separated, one row per stream in the order
Deg01..Deg07 then Pres01..Pres07, with columns for forest, cropland and urban
cover, land cover class diversity, the Human Footprint Index, elevation,
latitude and longitude.

## Running

```bash
BRACKEN=path/to/bracken
KO=path/to/Table_S4_KO_abundance.xlsx
LAND=LandUseCover1.csv

python3 scripts/01_genus_barplot.py     --bracken $BRACKEN
python3 scripts/02_alpha_diversity.py   --bracken $BRACKEN
python3 scripts/03_nmds_permanova.py    --bracken $BRACKEN --ko $KO
python3 scripts/04_guilds_redox.py      --bracken $BRACKEN
python3 scripts/05_marker_genes.py      --ko $KO --land $LAND
python3 scripts/06_dbrda_mantel.py      --bracken $BRACKEN --land $LAND
python3 scripts/07_spatial_confounding.py --bracken $BRACKEN --land $LAND
python3 scripts/08_guilds_gradient.py \
    --guilds guilds_redox.csv --tests guilds_redox_tests.csv --land $LAND
```

Script 8 reads the output of script 4, so the order matters there. The rest
are independent.

## Decisions worth knowing about

**Guilds, not genera.** With seven samples per group, the smallest p value a
Mann-Whitney test can return is 5.8e-4, reached only when the two groups
separate perfectly. Across 573 genera that floor corresponds to a corrected q
of 0.33, so no genus can reach significance in this design whatever its
effect size. The analysis proceeds at guild level, where eleven tests replace
573 and each category aggregates many genera. A null result at genus level is
a property of the sample size, not evidence that no genus differs.

**Guilds ordered by electron acceptor.** The order in which the guilds are
listed is the figure's argument: nitrate, then Fe(III), then sulfate, then
CO2 is the sequence in which a sediment loses electron acceptors as oxygen is
consumed. Sorting them by effect size or alphabetically would hide it.

**Proportions within domain.** Archaea are 0.15% to 0.86% of classified
reads here, so a bacterial and an archaeal guild expressed as percentages of
the whole community are not comparable. Each is a percentage of its own
domain. The two must not be added.

**Halophile and thermophile archaea excluded.** Neither has a plausible
population in a freshwater tropical stream, and both are read as
misclassification against a database whose archaeal representation is
dominated by cultured extremophiles. They took 9% to 66% of the archaeal
reads, a large and variable share that distorted every other proportion. The
exclusion changes the result and is stated rather than silent.

**A riparian index, not a land use index.** Land cover was extracted within a
30 m buffer of each stream, which at 30 m pixel resolution is about one pixel
either side of the channel. That measures the condition of the riparian strip,
not land use across the contributing catchment. The first principal component
of the four cover metrics is used as a single predictor; it should be read as
riparian cover, and the buffer does not delimit a catchment.

**The Human Footprint Index is not a predictor.** Land cover is among the
layers it is built from, so cropland and urban cover would be both predictors
and inputs. It is used to check that the locally built index describes the
same gradient - the two agree at rho = 0.776 - and as a descriptor of
anthropogenic pressure in the limnological analysis.

**One predictor, defined in advance.** The constrained ordination uses the
index alone. An earlier version fitted a model on the metrics that had passed
marginal tests and reported its R2 as 0.49; that is selection on the outcome
and inflates the number. The reported value is 0.189.

**Limnological variables are not corrected.** Correcting across six variables
and eleven guilds would remove any result at this sample size and would
suggest that absence of association had been established when it had merely
not been detected. The magnitudes carry the interpretation: correlations with
dissolved oxygen are below 0.22 in absolute value, which no correction would
make interpretable.

## Data availability

Raw reads: NCBI SRA, BioProject [accession].
Bracken reports, KO table and landscape data: [archive DOI].

## Citation

[citation]

## Licence

[licence]
