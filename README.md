This reproduces and furthers the data and scripts deposited on [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.z34tmpgbb) for the paper:

Chong, K.Y., Ng, W.Q., A.T.K. Yee, D.L. Yong. (2021) The community structure of bird assemblages on urban strangler figs. Biotropica 53: 255–265. [[website](https://onlinelibrary.wiley.com/doi/abs/10.1111/btp.12866)] [[pdf](https://www.dropbox.com/s/zlb8lu04gmfg5m3/chong%20et%20al%20in%20press%20the%20community%20structure%20of%20bird%20assemblages%20on%20urban%20strangler%20figs.pdf?dl=0)]

# Abstract

Figs have been regarded as keystone plant resources that support diverse tropical vertebrate frugivore communities. Planting or conserving large fig trees, such as stranglers, has therefore been proposed for enhancing urban biodiversity. We compared the diversity and community structure of bird assemblages on strangler figs with non-fig urban trees as well as between the fruiting and non-fruiting fig trees in an urban setting in Singapore. The total bird abundance across all the fig trees when in fruit was 4.5-fold higher than on non-fig trees and 3.5-fold higher than when the same fig trees were not fruiting, but only attracted two more species. On individual trees, after accounting for the presence of mistletoes, tree height, the area covered by buildings and road lane density, and distance to natural vegetation, mean diversity was not different between non-fig trees and fig trees when they were not in fruit. On the other hand, when fruiting, each fig tree on average had 1.4 more species, 3 more counts of non-native birds, and 1.6 more counts of insectivorous birds than when not fruiting. There was significant compositional turnover between non-fig trees and non-fruiting fig trees, while community dispersion was significantly lower among fig trees in fruit. Our results demonstrate that fig trees provide fruit and non-fruit resources for birds in an urban landscape but do not necessarily support a more diverse total bird assemblages than non-fig trees. Instead, bird communities on fruiting urban figs would be highly homogeneous and dominated by a few species.

# Methods

## Tree selection

Thirty pairs of fig [Ficus benjamina (Fb), F. elastica (Fe), F. microcarpa (Fm), F. religiosa (Fr), and F. superba (Fs)] and non-fig [Albizia saman (= Samanea saman; Ss), Peltophorum pterocarpum (Pp), Swietenia macrophylla (Sm), Pterocarpus indicus (Pi), and Khaya senegalensis (Ks)] focal trees were selected in this study. We first considered fig trees that were clearly in an urban landscape for our study, i.e., surrounded by built-up area and managed greenery but without natural vegetation within a 50-m radius from the tree’s trunk. The criteria for selecting fig focal trees included a minimum height (at least 10 m) which also ensures reproductive maturity. A non-fig focal tree was then identified near the potential fig focal tree. Non-fig focal trees were selected to match as closely as possible to its paired fig tree in terms of tree height as well as the urban landscape characteristics.

Seven fig trees did not display sufficient fruiting to qualify for a FIGF survey, while one F. elastica was continuously in the edible phase as defined above and did not have a FIGNF survey. See below for more information on FIGF and FIGNF in addition to NONFIG surveys.

## Bird surveys

Bird counts for each tree were conducted in every other 15-min interval from 7:30 am to 9:30 am, on days without strong wind or rain. All birds sighted (with the aid of an Opticron 8×42 binoculars) or heard within 20 m of the tree trunk by a single observer (Ng W.Q.) were identified and counted. Overflying birds were excluded. Birds seen or heard but could not be identified during the count were either photographed or the calls were recorded. These were later identified using appropriate references.Only one tree was surveyed at any one time. All bird counts were made from 19 August 2011 to 15 February 2012, which encompasses the peak period for visits by migratory birds in Singapore.

Bird counts were conducted once on each non-fig tree (NONFIG) and twice on each fig tree: once when the fig tree displayed edible (or ripe) fruit (FIGF) and once when it was not (FIGNF). FIGF bird counts were only conducted when edible figs on the tree were estimated to be at least 50% of the total crop produced. This was quantified by partitioning the tree crown into ten approximately equally-spaced wedges radiating from the trunk and estimating the percentage of ripe figs in each wedge, and then taking the average.

As there are some non-focal fig trees near our focal trees, no counts were made for the FIGNF and NONFIG survey types if such fig trees present within 200 m of their vicinity were producing edible fruits.

## Tree characteristics and urban landscape mapping

With an unobscured Google Earth image most recent to the survey date serving as base map, tree canopy cover (CANOPY; including both the focal tree crown as well as non-focal tree crowns), the area covered by buildings (BUILT) and by water bodies (WATER) were manually drawn at 1:1000 magnification. Most of the images were taken in the year of survey, i.e., 2011; however, older images from 2009 and 2010 were used if the latest image was obstructed by cloud cover.

This was then ground-truthed on the same day as the bird surveys, together with the data collection for other covariates: the heights of the focal tree (TREE_HT) and the tallest building within 50 m (BUILD_HT). In addition, the presence or absence of mistletoes, such as Dendrophthoe pentandra or Macrosolen cochinchinensis (Family: Loranthaceae), on the focal trees (MISTLETOE) was also recorded.

The distance of the focal tree to the nearest natural vegetation patch (NAT_DIST) was extracted from the vegetation map of Singapore (Yee et al. 2011 Gardens' Bulletin ).

# Usage notes

We have also provided the R script that produces the results in the paper. This R script uses each of the spreadsheets in the .XLSX file exported as .CSV files as inputs; the name of each spreadsheet is the name of the input .CSV file, e.g., "bird.csv".

The .ZIP file contains all the shapefiles and associated files from the urban landscape mapping. These also need to be unzipped into a folder and the folder path modified accordingly in the R script to be read and the relevant landscape variables extracted and calculated.

# Funding

Nature Society Singapore
