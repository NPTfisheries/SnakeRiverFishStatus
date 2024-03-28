
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SnakeRiverFishStatus <a href='https://github.com/NPTfisheries/SnakeRiverFishStatus'><img src='DFRM.png' align="right" width="110" /></a>

## Description

Welcome to the `SnakeRiverFishStatus` GitHub repository! The repo
contains data, scripts, functions, and workflow to estimate adult
Chinook salmon and steelhead escapement to Lower Granite Dam (LGR) as
well as to instream PIT-tag detection systems (IPTDS) throughout the
Snake River basin. Adult abundances are also provided for Snake River
Interior Columbia Technical Recovery Team (ICTRT) populations, where
available, and are the summation of one or more IPTDS sites within the
population.

Estimates of total adult escapement at LGR use the **ST**ate space
**A**dult **D**am **E**scapement **M**odel (STADEM) which is described
in [See et
al. (2021)](https://afspubs.onlinelibrary.wiley.com/doi/abs/10.1002/nafm.10649).
The STADEM R package is maintained
[here](https://github.com/KevinSee/STADEM). The movement probabilities
of upstream migrating adults to any “patch” (i.e., locations between or
past any IPTDS) throughout the basin is estimated using the **D**am
**A**dult **B**ranch **O**ccupancy **M**odel (DABOM) which is described
in [Waterhouse et
al. (2020)](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/eap.2202).
The DABOM R package is maintained
[here](https://github.com/KevinSee/DABOM). The `SnakeRiverFishStatus`
repo also widely leverages functionality in
[PITcleanr](https://github.com/KevinSee/PITcleanr) to map stream
networks throughout the basin and to wrangle large PIT tag observation
datasets from [PTAGIS](https://www.ptagis.org/). Posteriors from the
STADEM and DABOM models are combined to estimate abundances to IPTDS
sites and populations. Finally, sex, age, and size information from
adults sampled at LGR and detected within populations are used to parse
abundance estimates.

## Authors

The `SnakeRiverFishStatus` is a collaborative project, with the primary
contributors being:

- Mike Ackerman (Nez Perce Tribe - Department of Fisheries Resource
  Management)
- Ryan N. Kinzer (Nez Perce Tribe - Department of Fisheries Resource
  Management)
- Kevin See (Washington Department of Fish & Wildlife)

### Licenses

See the [LICENSE](LICENSE) file.
