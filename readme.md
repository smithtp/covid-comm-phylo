#  COVID community phylogenetics metrics

The COVID-19 pandemic provides an unprecedented opportunity to develop and test new metrics and methods for spotting disease outbreaks. Can classic metrics of ecological community structure (specifically, MPD and MNTD) tell us anything useful about how the COVID variants have spread through time and could these be a useful addition to traditional phylogenetic tools in studying outbreaks? 

## Repository Contents:

### Data 
All trees, meta data and dispersion metrics dataframes can be found in `data` 
folder. Manon Ragonnet-Cronin from the phylodynamics team investigating COVID
at Imperial generated the German, French, Italian and Spainish trees. The European tree was downloaded from NextStrain on 2021-04-26.

### Code

Scripts to run analysis and produce figures, in order:

 1. `phylo_shape.R` - Reads the phylogenies and creates community phylogenetic metrics. Script is separated by the different datasets, producing a `data/$_dispersion.csv` for each dataset (where `$` = Europe/Germany/Italy/Spain/France). All metadata is also collected in `data/metas.RData`.
 2. `data-wrangling.R` - creates additional variables for analyses. Reads the dispersion metrics and metadata files, produces `data/model_dfs.RData`.
 3. `phylo-results.R` - main analysis code. Reads `data/model_dfs.RData`, creates figures and runs GAM models.
 4. `maps.R` - plots spatial varition in dispersion metrics through time.
 5. `exploratory.R` - a boat-load of additional exploratory figures to investigate whats going on in the datasets. See `Jess_EcoCovid_Notes.pdf` for a more detailed explanation of what was being tested here.

### Phylogenetic dispersion metrics 
To re-calculate metrics run `code/phylo_shape.R`, some of these take a long time 
to run. From pez, '.obs.z' is the Standardised effect size; mpd.z = -NRI; mntd.z = -NTI. d.D = D
(NRI = "Net related index", NTI = "Nearest taxon index").

NRI & NTI:
 * Positive NRI & NTI = clustering = more species per genus than expected 
 * Negative NRI & NTI = overdispersed = less species per genus than expected 

D: 
 * 1 = Random with respect to phylo
 * 0 = Expected with Brownian motion
 * &gt; 1 = More overdispersion than expected by random
 * &lt; 0 = More clustered than expect under brownian


### Results

Both MPD and MNTD seem to respond strongly to the spread of the Alpha, and then Delta, showing a sharp increase in phylogenetic clustering with the Alpha wave, then a return towards the baseline with Delta (see figures).

Whats driving these differences in community structure through time?