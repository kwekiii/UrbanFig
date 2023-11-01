# load libraries
library(rgdal)
library(rgeos)
library(raster)
library(vegan)
library(lme4)
library(ape)
library(MuMIn)
library(MASS)
# library(ggplot2)
library(gridExtra)
library(betapart)
library(iNEXT)
library(labdsv)
library(tidyverse)
library(readxl)

# load landscape and tree data

sourceFilePath <- "data/UrbanFig_data.xlsx"

env <- read_excel(sourceFilePath, sheet = "env")
env$tree_type<-factor(env$tree_type, levels=levels(env$tree_type)[c(2,1)])

trees.ly<-readOGR(".//../GIS", "tree location")
canopy.ly<-readOGR(".//../GIS", "L126_canopycover")
canopy.ly<-gBuffer(canopy.ly, width = 0) #canopy layer contains self-intersecting geometries
roads.ly<-readOGR(".//../GIS", "L126_road")
others.ly<-readOGR(".//../GIS", "L126_built+waterbody")

## 126-m radius buffer

trees_buffer126<-gBuffer(trees.ly, width = 126, byid = TRUE)

canopy_buffer126<-intersect(canopy.ly, trees_buffer126)
env$canopy126<-gArea(canopy_buffer126, byid = TRUE)[match(env$tree_ID, canopy_buffer126$Name)]
env$canopy126<-env$canopy126/(pi*126^2)
env$canopy126<-replace(env$canopy126, is.na(env$canopy126), 0)

roads_buffer126<-intersect(roads.ly, trees_buffer126)
roads_buffer126@data$length<-gLength(roads_buffer126, byid = TRUE)
lengthxlanes<-with(roads_buffer126@data, lanes*length)
env$roads126<-tapply(lengthxlanes, roads_buffer126@data$Name.2, sum)/(pi*126^2)
env$roads126<-replace(env$roads126, is.na(env$roads126), 0)

others_buffer126<-intersect(others.ly, trees_buffer126)

water_buffer126<-others_buffer126[others_buffer126@data$type=="wb",]
water_buffer126@data$area<-gArea(water_buffer126, byid = TRUE)
env$water126<-with(water_buffer126@data, tapply(area, Name.2, sum))/(pi*126^2)
env$water126<-replace(env$water126, is.na(env$water126), 0)

built_buffer126<-others_buffer126[others_buffer126@data$type=="b",]
built_buffer126@data$area<-gArea(built_buffer126, byid = TRUE)
env$built126<-with(built_buffer126@data, tapply(area, Name.2, sum))/(pi*126^2)
env$built126<-replace(env$built126, is.na(env$built126), 0)

## 50-m radius buffer

trees_buffer50<-gBuffer(trees.ly, width = 50, byid = TRUE)

canopy_buffer50<-intersect(canopy.ly, trees_buffer50)
env$canopy50<-gArea(canopy_buffer50, byid = TRUE)[match(env$tree_ID, canopy_buffer50$Name)]
env$canopy50<-env$canopy50/(pi*50^2)
env$canopy50<-replace(env$canopy50, is.na(env$canopy50), 0)

roads_buffer50<-intersect(roads.ly, trees_buffer50)
roads_buffer50@data$length<-gLength(roads_buffer50, byid = TRUE)
lengthxlanes<-with(roads_buffer50@data, lanes*length)
env$roads50<-tapply(lengthxlanes, roads_buffer50@data$Name.2, sum)/(pi*50^2)
env$roads50<-replace(env$roads50, is.na(env$roads50), 0)

others_buffer50<-intersect(others.ly, trees_buffer50)

water_buffer50<-others_buffer50[others_buffer50@data$type=="wb",]
water_buffer50@data$area<-gArea(water_buffer50, byid = TRUE)
env$water50<-with(water_buffer50@data, tapply(area, Name.2, sum))/(pi*50^2)
env$water50<-replace(env$water50, is.na(env$water50), 0)

built_buffer50<-others_buffer50[others_buffer50@data$type=="b",]
built_buffer50@data$area<-gArea(built_buffer50, byid = TRUE)
env$built50<-with(built_buffer50@data, tapply(area, Name.2, sum))/(pi*50^2)
env$built50<-replace(env$built50, is.na(env$built50), 0)

env.unscaled<-env
env<-env[,-c(6,10:13)]
env[,c(6:16)]<-scale(env.unscaled[,c(7:9,14:21)])

# load bird survey data

bird <- read_excel(sourceFilePath, sheet = "bird")
status <- read_excel(sourceFilePath, sheet = "status")

## nonfig vs. fignf
nfvfnf.com<-xtabs(count~tree_ID+bird_sp, data=
	aggregate(count~tree_ID+bird_sp, FUN=max, data=bird[bird$survey_type!="figf",])
	)
nfvfnf.com<-as.data.frame.matrix(nfvfnf.com)

nfvfnf.D0<-specnumber(nfvfnf.com)
nfvfnf.D1<-exp(diversity(nfvfnf.com, index="shannon"))
nfvfnf.D2<-diversity(nfvfnf.com, index="invsimpson")

nfvfnf.e.com<-nfvfnf.com[,status$status=="e"]
nfvfnf.eabund<-apply(nfvfnf.e.com,1,sum)

nfvfnf.I.com<-nfvfnf.com[,status$I==1&status$F==0]
nfvfnf.Iabund<-apply(nfvfnf.I.com,1,sum)

## fignf vs. figf
bird$index2<-with(bird, paste(survey_type, tree_ID))
fnfvff.com<-xtabs(count~index2+bird_sp, data=
	aggregate(count~index2+bird_sp, FUN=max, data=bird[bird$survey_type!="nonfig",])
	)
fnfvff.com<-as.data.frame.matrix(fnfvff.com)

fnfvff.D0<-specnumber(fnfvff.com)
fnfvff.D1<-exp(diversity(fnfvff.com, index="shannon"))
fnfvff.D2<-diversity(fnfvff.com, index="invsimpson")

fnfvff.e.com<-fnfvff.com[,status$status=="e"]
fnfvff.eabund<-apply(fnfvff.e.com,1,sum)

fnfvff.I.com<-fnfvff.com[,status$I==1&status$F==0]
fnfvff.Iabund<-apply(fnfvff.I.com,1,sum)

figging<-sapply(strsplit(rownames(fnfvff.com), " "), "[")[1,]
figging<-factor(figging, levels=c("fignf","figf"))

tree_ID2<-sapply(strsplit(rownames(fnfvff.com), " "), "[")[2,]
tree_sp2<-env$tree_sp[match(tree_ID2, env$tree_ID)]

# GLMMs

## nonfig vs. fignf

### Species richness

nfvfnf.D0.mod.null1<-glmer(nfvfnf.D0~1+(1|tree_sp)+(1|site),
                           data=env, subset=env$tree_ID!="D2Fe", family=poisson, na.action="na.fail")
nfvfnf.D0.mod.null2<-glmer(nfvfnf.D0~1+(1|tree_sp),
                           data=env, subset=env$tree_ID!="D2Fe", family=poisson, na.action="na.fail")
nfvfnf.D0.mod.null3<-glmer(nfvfnf.D0~1+(1|site),
                           data=env, subset=env$tree_ID!="D2Fe", family=poisson, na.action="na.fail")
nfvfnf.D0.mod.null4<-glm(nfvfnf.D0~1,
                           data=env, subset=env$tree_ID!="D2Fe", family=poisson, na.action="na.fail")

nfvfnf.D0.mod.full<-glm(nfvfnf.D0~tree_type+nat_dist+tree_ht+build_ht+mistletoe+built50+built126+canopy50+canopy126+water50+water126+roads50+roads126,
                          data=env, subset=env$tree_ID!="D2Fe", family=poisson, na.action="na.fail")
nfvfnf.D0.mod.tab<-dredge(nfvfnf.D0.mod.full,
       subset= !(built50 & built126) & !(canopy50 & canopy126) & !(water50 & water126) & !(roads50 & roads126),
       m.lim=c(NA, 6), extra = "R^2")

nfvfnf.D0.mod.tab.top15<-data.frame(head(nfvfnf.D0.mod.tab,15))
colnames(nfvfnf.D0.mod.tab.top15)[1]<-"(Intercept)"
colnames(nfvfnf.D0.mod.tab.top15)[12]<-"tree_typefig"
nfvfnf.D0.mod.tab.top15$tree_typefig<-as.numeric(nfvfnf.D0.mod.tab.top15$tree_typefig)

nfvfnf.D0.mod.tab.top15.se<-matrix(NA, nrow = nrow(nfvfnf.D0.mod.tab.top15), ncol = 14)
colnames(nfvfnf.D0.mod.tab.top15.se)<-colnames(nfvfnf.D0.mod.tab.top15)[1:14]

for (i in 1:nrow(nfvfnf.D0.mod.tab.top15)) {
  coef.tab<-coefTable(nfvfnf.D0.mod.tab[i])[[1]]
  nfvfnf.D0.mod.tab.top15[i,match(rownames(coef.tab), colnames(nfvfnf.D0.mod.tab.top15))]<-coef.tab[,1]
  nfvfnf.D0.mod.tab.top15.se[i,match(rownames(coef.tab), colnames(nfvfnf.D0.mod.tab.top15.se))]<-coef.tab[,2]
}

write.csv(
  cbind(nfvfnf.D0.mod.tab.top15, nfvfnf.D0.mod.tab.top15.se)[,c(as.vector(rbind(1:14,21:34)),15:20)]
  , "nfvfnf_D0_modtab.csv", quote = FALSE, row.names = FALSE
)

nfvfnf.D0.mod.top1<-get.models(nfvfnf.D0.mod.tab, subset = 1)
deviance(nfvfnf.D0.mod.top1[[1]])/df.residual(nfvfnf.D0.mod.top1[[1]]) # no overdispersion

### First-order diversity

nfvfnf.D1.mod.null1<-lmer(log(nfvfnf.D1)~1+(1|tree_sp)+(1|site),
                      data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D1.mod.null2<-lmer(log(nfvfnf.D1)~1+(1|tree_sp),
                      data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D1.mod.null3<-lmer(log(nfvfnf.D1)~1+(1|site),
                      data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D1.mod.null4<-lm(log(nfvfnf.D1)~1,
                    data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")

nfvfnf.D1.mod.full<-lm(log(nfvfnf.D1)~tree_type+nat_dist+tree_ht+build_ht+mistletoe+built50+built126+canopy50+canopy126+water50+water126+roads50+roads126,
                          data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D1.mod.tab<-dredge(nfvfnf.D1.mod.full,
                          subset= !(built50 & built126) & !(canopy50 & canopy126) & !(water50 & water126) & !(roads50 & roads126),
                          m.lim = c(NA, 6), extra = "R^2")

nfvfnf.D1.mod.tab.top25<-data.frame(head(nfvfnf.D1.mod.tab,25))
colnames(nfvfnf.D1.mod.tab.top25)[1]<-"(Intercept)"
colnames(nfvfnf.D1.mod.tab.top25)[12]<-"tree_typefig"
nfvfnf.D1.mod.tab.top25$tree_typefig<-as.numeric(nfvfnf.D1.mod.tab.top25$tree_typefig)

nfvfnf.D1.mod.tab.top25.se<-matrix(NA, nrow = nrow(nfvfnf.D1.mod.tab.top25), ncol = 14)
colnames(nfvfnf.D1.mod.tab.top25.se)<-colnames(nfvfnf.D1.mod.tab.top25)[1:14]

for (i in 1:nrow(nfvfnf.D1.mod.tab.top25)) {
  coef.tab<-coefTable(nfvfnf.D1.mod.tab[i])[[1]]
  nfvfnf.D1.mod.tab.top25[i,match(rownames(coef.tab), colnames(nfvfnf.D1.mod.tab.top25))]<-coef.tab[,1]
  nfvfnf.D1.mod.tab.top25.se[i,match(rownames(coef.tab), colnames(nfvfnf.D1.mod.tab.top25.se))]<-coef.tab[,2]
}

write.csv(
  cbind(nfvfnf.D1.mod.tab.top25, nfvfnf.D1.mod.tab.top25.se)[,c(as.vector(rbind(1:14,21:34)),15:20)]
  , "nfvfnf_D1_modtab.csv", quote = FALSE, row.names = FALSE
)

nfvfnf.D1.mod.top1<-get.models(nfvfnf.D1.mod.tab, subset = 1)

### Second-order diversity

nfvfnf.D2.mod.null1<-lmer(log(nfvfnf.D2)~1+(1|tree_sp)+(1|site),
                     data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D2.mod.null2<-lmer(log(nfvfnf.D2)~1+(1|tree_sp),
                     data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D2.mod.null3<-lmer(log(nfvfnf.D2)~1+(1|site),
                     data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D2.mod.null4<-lm(log(nfvfnf.D2)~1,
                   data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")

nfvfnf.D2.mod.full<-lm(log(nfvfnf.D2)~tree_type+nat_dist+tree_ht+build_ht+mistletoe+built50+built126+canopy50+canopy126+water50+water126+roads50+roads126,
                          data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.D2.mod.tab<-dredge(nfvfnf.D2.mod.full,
                          subset= !(built50 & built126) & !(canopy50 & canopy126) & !(water50 & water126) & !(roads50 & roads126),
                          m.lim = c(NA, 6), extra = "R^2")

nfvfnf.D2.mod.tab.top30<-data.frame(head(nfvfnf.D2.mod.tab,30))
colnames(nfvfnf.D2.mod.tab.top30)[1]<-"(Intercept)"
colnames(nfvfnf.D2.mod.tab.top30)[12]<-"tree_typefig"
nfvfnf.D2.mod.tab.top30$tree_typefig<-as.numeric(nfvfnf.D2.mod.tab.top30$tree_typefig)

nfvfnf.D2.mod.tab.top30.se<-matrix(NA, nrow = nrow(nfvfnf.D2.mod.tab.top30), ncol = 14)
colnames(nfvfnf.D2.mod.tab.top30.se)<-colnames(nfvfnf.D2.mod.tab.top30)[1:14]

for (i in 1:nrow(nfvfnf.D2.mod.tab.top30)) {
  coef.tab<-coefTable(nfvfnf.D2.mod.tab[i])[[1]]
  nfvfnf.D2.mod.tab.top30[i,match(rownames(coef.tab), colnames(nfvfnf.D2.mod.tab.top30))]<-coef.tab[,1]
  nfvfnf.D2.mod.tab.top30.se[i,match(rownames(coef.tab), colnames(nfvfnf.D2.mod.tab.top30.se))]<-coef.tab[,2]
}

write.csv(
  cbind(nfvfnf.D2.mod.tab.top30, nfvfnf.D2.mod.tab.top30.se)[,c(as.vector(rbind(1:14,21:34)),15:20)]
  , "nfvfnf_D2_modtab.csv", quote = FALSE, row.names = FALSE
)

nfvfnf.D2.mod.top1<-get.models(nfvfnf.D2.mod.tab, subset = 1)

### Exotic abundance

nfvfnf.eabund.mod.null1<-glmer.nb(nfvfnf.eabund~1+(1|tree_sp)+(1|site),
                               data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.eabund.mod.null2<-glmer.nb(nfvfnf.eabund~1+(1|tree_sp),
                                  data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.eabund.mod.null3<-glmer.nb(nfvfnf.eabund~1+(1|site),
                                  data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.eabund.mod.null4<-glm.nb(nfvfnf.eabund~1,
                                  data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")

nfvfnf.eabund.mod.full<-glm.nb(nfvfnf.eabund~tree_type+nat_dist+tree_ht+build_ht+mistletoe+built50+built126+canopy50+canopy126+water50+water126+roads50+roads126,
                       data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.eabund.mod.tab<-dredge(nfvfnf.eabund.mod.full,
                          subset= !(built50 & built126) & !(canopy50 & canopy126) & !(water50 & water126) & !(roads50 & roads126),
                          m.lim = c(NA, 6), extra = "R^2")

nfvfnf.eabund.mod.tab.top10<-data.frame(head(nfvfnf.eabund.mod.tab,10))
colnames(nfvfnf.eabund.mod.tab.top10)[1]<-"(Intercept)"
colnames(nfvfnf.eabund.mod.tab.top10)[12]<-"tree_typefig"
nfvfnf.eabund.mod.tab.top10$tree_typefig<-as.numeric(nfvfnf.eabund.mod.tab.top10$tree_typefig)

nfvfnf.eabund.mod.tab.top10.se<-matrix(NA, nrow = nrow(nfvfnf.eabund.mod.tab.top10), ncol = 14)
colnames(nfvfnf.eabund.mod.tab.top10.se)<-colnames(nfvfnf.eabund.mod.tab.top10)[1:14]

for (i in 1:nrow(nfvfnf.eabund.mod.tab.top10)) {
  coef.tab<-coefTable(nfvfnf.eabund.mod.tab[i])[[1]]
  nfvfnf.eabund.mod.tab.top10[i,match(rownames(coef.tab), colnames(nfvfnf.eabund.mod.tab.top10))]<-coef.tab[,1]
  nfvfnf.eabund.mod.tab.top10.se[i,match(rownames(coef.tab), colnames(nfvfnf.eabund.mod.tab.top10.se))]<-coef.tab[,2]
}

write.csv(
  cbind(nfvfnf.eabund.mod.tab.top10, nfvfnf.eabund.mod.tab.top10.se)[,c(as.vector(rbind(1:14,21:34)),15:20)]
  , "nfvfnf_eabund_modtab.csv", quote = FALSE, row.names = FALSE
)

nfvfnf.eabund.mod.top1<-get.models(nfvfnf.eabund.mod.tab, subset = 1)

### Insectivore abundance

nfvfnf.Iabund.mod.null1<-glmer.nb(nfvfnf.Iabund~1+(1|tree_sp)+(1|site), control = glmerControl(optimizer="bobyqa"),
                                  data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.Iabund.mod.null2<-glmer.nb(nfvfnf.Iabund~1+(1|tree_sp),
                                  data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.Iabund.mod.null3<-glmer.nb(nfvfnf.Iabund~1+(1|site),
                                  data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.Iabund.mod.null4<-glm.nb(nfvfnf.Iabund~1,
                                data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")

nfvfnf.Iabund.mod.full<-glm.nb(nfvfnf.Iabund~tree_type+nat_dist+tree_ht+build_ht+mistletoe+built50+built126+canopy50+canopy126+water50+water126+roads50+roads126,
                               data=env, subset=env$tree_ID!="D2Fe", na.action="na.fail")
nfvfnf.Iabund.mod.tab<-dredge(nfvfnf.Iabund.mod.full,
                              subset= !(built50 & built126) & !(canopy50 & canopy126) & !(water50 & water126) & !(roads50 & roads126),
                              m.lim = c(NA, 6), extra = "R^2")

nfvfnf.Iabund.mod.tab.top20<-data.frame(head(nfvfnf.Iabund.mod.tab,20))
colnames(nfvfnf.Iabund.mod.tab.top20)[1]<-"(Intercept)"
colnames(nfvfnf.Iabund.mod.tab.top20)[12]<-"tree_typefig"
nfvfnf.Iabund.mod.tab.top20$tree_typefig<-as.numeric(nfvfnf.Iabund.mod.tab.top20$tree_typefig)

nfvfnf.Iabund.mod.tab.top20.se<-matrix(NA, nrow = nrow(nfvfnf.Iabund.mod.tab.top20), ncol = 14)
colnames(nfvfnf.Iabund.mod.tab.top20.se)<-colnames(nfvfnf.Iabund.mod.tab.top20)[1:14]

for (i in 1:nrow(nfvfnf.Iabund.mod.tab.top20)) {
  coef.tab<-coefTable(nfvfnf.Iabund.mod.tab[i])[[1]]
  nfvfnf.Iabund.mod.tab.top20[i,match(rownames(coef.tab), colnames(nfvfnf.Iabund.mod.tab.top20))]<-coef.tab[,1]
  nfvfnf.Iabund.mod.tab.top20.se[i,match(rownames(coef.tab), colnames(nfvfnf.Iabund.mod.tab.top20.se))]<-coef.tab[,2]
}

write.csv(
  cbind(nfvfnf.Iabund.mod.tab.top20, nfvfnf.Iabund.mod.tab.top20.se)[,c(as.vector(rbind(1:14,21:34)),15:20)]
  , "nfvfnf_Iabund_modtab.csv", quote = FALSE, row.names = FALSE
)

nfvfnf.Iabund.mod.top1<-get.models(nfvfnf.Iabund.mod.tab, subset = 1)

## fignf vs. figf

### Species richness

fnfvff.D0.mod.null1<-glmer(fnfvff.D0~1+(1|tree_ID2)+(1|tree_sp2), family=poisson)
fnfvff.D0.mod.null2<-glmer(fnfvff.D0~1+(1|tree_ID2), family=poisson)

fnfvff.D0.mod<-glmer(fnfvff.D0~figging+(1|tree_ID2), family=poisson)

fnfvff.D0.mod.tab<-model.sel(list(fnfvff.D0.mod.null2, fnfvff.D0.mod), extra = list(R2 = function(x) r.squaredGLMM(x, fnfvff.D0.mod.null2)["trigamma",]))
fnfvff.D0.mod.tab1<-data.frame(fnfvff.D0.mod.tab)
fnfvff.D0.mod.tab1$figging<-as.numeric(fnfvff.D0.mod.tab1$figging)
colnames(fnfvff.D0.mod.tab1)[1:2]<-c("(Intercept)", "figgingfigf")

fnfvff.D0.mod.tab.se<-data.frame(matrix(NA, nrow=2, ncol=2))
colnames(fnfvff.D0.mod.tab.se) <- c("(Intercept)", "figgingfigf")

for (i in 1:2) {
  coef.tab<-coefTable(fnfvff.D0.mod.tab)[[i]]
  fnfvff.D0.mod.tab1[i, match(rownames(coef.tab), colnames(fnfvff.D0.mod.tab1))]<-coef.tab[,1]
  fnfvff.D0.mod.tab.se[i, match(rownames(coef.tab), colnames(fnfvff.D0.mod.tab.se))]<-coef.tab[,2]
}

fnfvff.D0.mod.tab<-cbind(fnfvff.D0.mod.tab1, fnfvff.D0.mod.tab.se)[,c(1,10,2,11,3:9)]

deviance(fnfvff.D0.mod)/df.residual(fnfvff.D0.mod) #under-dispersed

### First order diversity

fnfvff.D1.mod.null1<-lmer(log(fnfvff.D1)~1+(1|tree_ID2)+(1|tree_sp2))
fnfvff.D1.mod.null2<-lmer(log(fnfvff.D1)~1+(1|tree_ID2))

fnfvff.D1.mod.null2a<-lmer(log(fnfvff.D1)~1+(1|tree_ID2), REML = FALSE)
fnfvff.D1.mod<-lmer(log(fnfvff.D1)~figging+(1|tree_ID2), REML = FALSE)

fnfvff.D1.mod.tab<-model.sel(list(fnfvff.D1.mod.null2a, fnfvff.D1.mod), extra = r.squaredGLMM)
fnfvff.D1.mod.tab1<-data.frame(fnfvff.D1.mod.tab)
fnfvff.D1.mod.tab1$figging<-as.numeric(fnfvff.D1.mod.tab1$figging)
colnames(fnfvff.D1.mod.tab1)[1:4]<-c("(Intercept)", "figgingfigf", "R2.R2m", "R2.R2c")

fnfvff.D1.mod.tab.se<-data.frame(matrix(NA, nrow=2, ncol=2))
colnames(fnfvff.D1.mod.tab.se) <- c("(Intercept)", "figgingfigf")

for (i in 1:2) {
  coef.tab<-coefTable(fnfvff.D1.mod.tab)[[i]]
  fnfvff.D1.mod.tab1[i, match(rownames(coef.tab), colnames(fnfvff.D1.mod.tab1))]<-coef.tab[,1]
  fnfvff.D1.mod.tab.se[i, match(rownames(coef.tab), colnames(fnfvff.D1.mod.tab.se))]<-coef.tab[,2]
}

fnfvff.D1.mod.tab<-cbind(fnfvff.D1.mod.tab1, fnfvff.D1.mod.tab.se)[,c(1,10,2,11,3:9)]

### Second order diversity

fnfvff.D2.mod.null1<-lmer(log(fnfvff.D2)~1+(1|tree_ID2)+(1|tree_sp2))
fnfvff.D2.mod.null2<-lmer(log(fnfvff.D2)~1+(1|tree_ID2))

fnfvff.D2.mod.null2a<-lmer(log(fnfvff.D2)~1+(1|tree_ID2), REML = FALSE)
fnfvff.D2.mod<-lmer(log(fnfvff.D2)~figging+(1|tree_ID2), REML = FALSE)

fnfvff.D2.mod.tab<-model.sel(list(fnfvff.D2.mod.null2a, fnfvff.D2.mod), extra = r.squaredGLMM)
fnfvff.D2.mod.tab1<-data.frame(fnfvff.D2.mod.tab)
fnfvff.D2.mod.tab1$figging<-as.numeric(fnfvff.D2.mod.tab1$figging)
colnames(fnfvff.D2.mod.tab1)[1:4]<-c("(Intercept)", "figgingfigf", "R2.R2m", "R2.R2c")

fnfvff.D2.mod.tab.se<-data.frame(matrix(NA, nrow=2, ncol=2))
colnames(fnfvff.D2.mod.tab.se) <- c("(Intercept)", "figgingfigf")

for (i in 1:2) {
  coef.tab<-coefTable(fnfvff.D2.mod.tab)[[i]]
  fnfvff.D2.mod.tab1[i, match(rownames(coef.tab), colnames(fnfvff.D2.mod.tab1))]<-coef.tab[,1]
  fnfvff.D2.mod.tab.se[i, match(rownames(coef.tab), colnames(fnfvff.D2.mod.tab.se))]<-coef.tab[,2]
}

fnfvff.D2.mod.tab<-cbind(fnfvff.D2.mod.tab1, fnfvff.D2.mod.tab.se)[,c(1,10,2,11,3:9)]

### Exotic abundance

fnfvff.eabund.mod.null1<-glmer.nb(fnfvff.eabund~1+(1|tree_ID2)+(1|tree_sp2), control = glmerControl(optimizer = "bobyqa"))
fnfvff.eabund.mod.null2<-glmer.nb(fnfvff.eabund~1+(1|tree_ID2))

fnfvff.eabund.mod<-glmer.nb(fnfvff.eabund~figging+(1|tree_ID2))

fnfvff.eabund.mod.tab<-model.sel(list(fnfvff.eabund.mod.null2, fnfvff.eabund.mod), extra = list(R2 = function(x) r.squaredGLMM(x, fnfvff.eabund.mod.null2)["trigamma",]))
fnfvff.eabund.mod.tab1<-data.frame(fnfvff.eabund.mod.tab)
fnfvff.eabund.mod.tab1$figging<-as.numeric(fnfvff.eabund.mod.tab1$figging)
colnames(fnfvff.eabund.mod.tab1)[1:2]<-c("(Intercept)", "figgingfigf")

fnfvff.eabund.mod.tab.se<-data.frame(matrix(NA, nrow=2, ncol=2))
colnames(fnfvff.eabund.mod.tab.se) <- c("(Intercept)", "figgingfigf")

for (i in 1:2) {
  coef.tab<-coefTable(fnfvff.eabund.mod.tab)[[i]]
  fnfvff.eabund.mod.tab1[i, match(rownames(coef.tab), colnames(fnfvff.eabund.mod.tab1))]<-coef.tab[,1]
  fnfvff.eabund.mod.tab.se[i, match(rownames(coef.tab), colnames(fnfvff.eabund.mod.tab.se))]<-coef.tab[,2]
}

fnfvff.eabund.mod.tab<-cbind(fnfvff.eabund.mod.tab1, fnfvff.eabund.mod.tab.se)[,c(1,11,2,12,3:4,6:10)]

### Insectivore abundance

fnfvff.Iabund.mod.null1<-glmer(fnfvff.Iabund~(1|tree_ID2)+(1|tree_sp2), family = poisson)
fnfvff.Iabund.mod.null2<-glmer(fnfvff.Iabund~(1|tree_ID2), family=poisson)

fnfvff.Iabund.mod<-glmer(fnfvff.Iabund~figging+(1|tree_ID2)+(1|tree_sp2), family = poisson)

deviance(fnfvff.Iabund.mod)/df.residual(fnfvff.Iabund.mod)

(fnfvff.Iabund.mod.tab<-model.sel(list(fnfvff.Iabund.mod.null1, fnfvff.Iabund.mod), extra = list(R2 = function(x) r.squaredGLMM(x, fnfvff.Iabund.mod.null1)["trigamma",])))
fnfvff.Iabund.mod.tab1<-data.frame(fnfvff.Iabund.mod.tab)
fnfvff.Iabund.mod.tab1$figging<-as.numeric(fnfvff.Iabund.mod.tab1$figging)
colnames(fnfvff.Iabund.mod.tab1)[1:2]<-c("(Intercept)", "figgingfigf")

fnfvff.Iabund.mod.tab.se<-data.frame(matrix(NA, nrow=2, ncol=2))
colnames(fnfvff.Iabund.mod.tab.se) <- c("(Intercept)", "figgingfigf")

for (i in 1:2) {
  coef.tab<-coefTable(fnfvff.Iabund.mod.tab)[[i]]
  fnfvff.Iabund.mod.tab1[i, match(rownames(coef.tab), colnames(fnfvff.Iabund.mod.tab1))]<-coef.tab[,1]
  fnfvff.Iabund.mod.tab.se[i, match(rownames(coef.tab), colnames(fnfvff.Iabund.mod.tab.se))]<-coef.tab[,2]
}

fnfvff.Iabund.mod.tab<-cbind(fnfvff.Iabund.mod.tab1, fnfvff.Iabund.mod.tab.se)[,c(1,10,2,11,3:9)]

### consolidated table

write.csv(
  rbind(fnfvff.D0.mod.tab, fnfvff.D1.mod.tab, fnfvff.D2.mod.tab, fnfvff.eabund.mod.tab, fnfvff.Iabund.mod.tab)
  , "fnfvff_modtab.csv", quote = FALSE
)

## other tables

nfvfnf.nulltab<-data.frame(
  D0=AICc(nfvfnf.D0.mod.null1, nfvfnf.D0.mod.null2, nfvfnf.D0.mod.null3, nfvfnf.D0.mod.null4)$AICc,
  D1=AICc(nfvfnf.D1.mod.null1, nfvfnf.D1.mod.null2, nfvfnf.D1.mod.null3, nfvfnf.D1.mod.null4)$AICc,
  D2=AICc(nfvfnf.D2.mod.null1, nfvfnf.D2.mod.null2, nfvfnf.D2.mod.null3, nfvfnf.D2.mod.null4)$AICc,
  eabund=AICc(nfvfnf.eabund.mod.null1, nfvfnf.eabund.mod.null2, nfvfnf.eabund.mod.null3, nfvfnf.eabund.mod.null4)$AICc,
  Iabund=AICc(nfvfnf.Iabund.mod.null1, nfvfnf.Iabund.mod.null2, nfvfnf.Iabund.mod.null3, nfvfnf.Iabund.mod.null4)$AICc,
  row.names = c("Species + site", "Species only", "Site only", "none"))

fnfvff.nulltab<-data.frame(
  rbind(
    AICc(fnfvff.D0.mod.null1, fnfvff.D1.mod.null1, fnfvff.D2.mod.null1, fnfvff.eabund.mod.null1, fnfvff.Iabund.mod.null1)$AICc,
    AICc(fnfvff.D0.mod.null2, fnfvff.D1.mod.null2, fnfvff.D2.mod.null2, fnfvff.eabund.mod.null2, fnfvff.Iabund.mod.null2)$AICc
  )
  , row.names = c("Species + tree", "Tree only")
)
colnames(fnfvff.nulltab)<-colnames(nfvfnf.nulltab)

write.csv(
 rbind(nfvfnf.nulltab, fnfvff.nulltab)
  , "AIC_nullmods.csv", quote = FALSE
)

## figures

boot.ci<-function(x, nperm=999, alpha=0.05, which) {
  xperm<-rep(NA, nperm)
  for(i in 1:nperm) {
    sample.x<-sample(x, length(x), replace=TRUE)
    xperm[i]<-mean(sample.x)
  }
  xperm<-sort(xperm)
  if(which=="lower") xperm[alpha/2*(nperm+1)] else
    if(which=="upper") xperm[(1-alpha/2)*(nperm+1)]
}

y.plot<-function(y.value, y.lab, tag, breaks) {
  df<-data.frame(
    type=c(as.character(env$tree_type[env$tree_ID!="D2Fe"]), as.character(figging[figging=="figf"])),
    tree_ID=c(as.character(env$tree_ID[env$tree_ID!="D2Fe"]), tree_ID2[figging=="figf"]),
    y.value
  )
  levels(df$type)<-c("FIGNF", "FIGF", "NONFIG")
  df$type<-factor(df$type, levels=c("NONFIG", "FIGNF", "FIGF"))
  
  df.tab<-aggregate(y.value~type, data=df, FUN=mean)
  names(df.tab)[2]<-"mean"
  df.tab$lower<-aggregate(y.value~type, data=df, FUN=boot.ci, which="lower")$y.value
  df.tab$upper<-aggregate(y.value~type, data=df, FUN=boot.ci, which="upper")$y.value
    
  g1<-ggplot(data=df.tab, aes(x=type, y=mean))+
    geom_bar(width = 0.7, stat="identity", alpha = c(1,1,0), fill=c("lightblue", "green", "pink"))
  
  g1wjitter<-g1+geom_jitter(data=df, mapping = aes(x=type, y=y.value), inherit.aes = FALSE, position = position_jitter(0.1), alpha=0.5, col="gray")
  g1wjitter.tab<-layer_data(g1wjitter, i=2)
  
  for_lines<-cbind(g1wjitter.tab, tree_ID=df$tree_ID)
  fig.both<-with(df, table(tree_ID,type))[,2:3]
  for_lines<-for_lines[df$tree_ID %in% rownames(fig.both[apply(fig.both, 1, sum)==2,]),]
  
  for_lines.tab<-as.data.frame.matrix(cbind(xtabs(x~droplevels(tree_ID)+group, data=for_lines), xtabs(y~droplevels(tree_ID)+group, data=for_lines)))
  colnames(for_lines.tab)<-c("x2","x3","y2","y3")
  
  g1 +
    geom_segment(data=for_lines.tab, aes(x = x2+0.1,
                                         y = y2,
                                         xend = x3-0.1,
                                         yend = y3), color="pink", size=1) +
    geom_point(data=g1wjitter.tab, mapping = aes(x=x, y=y), inherit.aes = FALSE, alpha=0.5, col=col.tree_type2, size=2, pch=pch.tree_type2)+
    labs(y = y.lab, x = "", tag = tag) + coord_cartesian(clip = "off") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(face = "italic"), 
          plot.tag = element_text(size = 10), 
          legend.position = "none", 
          plot.margin = margin(1.5,1.5,0.5,1.5))  
}

grid.arrange(
  y.plot(y.value = c(nfvfnf.D0[env$tree_ID!="D2Fe"], fnfvff.D0[figging=="figf"]),
       y.lab = "D0", tag = "(a)", breaks = c(-2, 0, 2, 4, 6)),
  y.plot(y.value = c(nfvfnf.D1[env$tree_ID!="D2Fe"], fnfvff.D1[figging=="figf"]),
         y.lab = "D1", tag = "(b)", breaks = c(-4, -2, 0, 2, 4)),
  y.plot(y.value = c(nfvfnf.D2[env$tree_ID!="D2Fe"], fnfvff.D2[figging=="figf"]),
         y.lab = "D2", tag = "(c)", breaks = c(-4, -2, 0, 2, 4)),
  y.plot(y.value = c(nfvfnf.eabund[env$tree_ID!="D2Fe"], fnfvff.eabund[figging=="figf"]),
         y.lab = "Exotic abundance", tag = "(d)", breaks = c(0, 10, 20, 30)),
  y.plot(y.value = c(nfvfnf.Iabund[env$tree_ID!="D2Fe"], fnfvff.Iabund[figging=="figf"]),
         y.lab = "Insectivore abundance", tag = "(e)", breaks = c(-2, 0, 2, 4)),
  ncol=3
)

# Beta diversity and NMDS

## nonfig vs. fignf

nfvfnf.bray<-bray.part(nfvfnf.com[env$tree_ID!="D2Fe",])

nfvfnf.com2<-ifelse(nfvfnf.com>0, 1, 0)
nfvfnf.beta<-beta.pair(nfvfnf.com2[env$tree_ID!="D2Fe",])

### with abundances

adonis(nfvfnf.bray$bray~mistletoe+nat_dist+tree_ht+build_ht, data=env[env$tree_ID!="D2Fe",], by = "margin")
#nat_dist N.S.

adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+canopy50, data=env[env$tree_ID!="D2Fe",], by = "margin") #marginally significant
adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+canopy126, data=env[env$tree_ID!="D2Fe",], by = "margin")

adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+roads50, data=env[env$tree_ID!="D2Fe",], by = "margin")
adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+roads126, data=env[env$tree_ID!="D2Fe",], by = "margin") #N.S.

adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+water50, data=env[env$tree_ID!="D2Fe",], by = "margin")
adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+water126, data=env[env$tree_ID!="D2Fe",], by = "margin") #marginally significant

adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+built50, data=env[env$tree_ID!="D2Fe",], by = "margin") #N.S.
adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+built126, data=env[env$tree_ID!="D2Fe",], by = "margin") #N.S.

adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+canopy126+roads50+water50, data=env[env$tree_ID!="D2Fe",], by = "margin") #water50 N.S.
adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+canopy126, data=env[env$tree_ID!="D2Fe",], by = "margin") #roads50 N.S.

### presence-absence only

adonis(nfvfnf.beta$beta.sor~mistletoe+nat_dist+tree_ht+build_ht, data=env[env$tree_ID!="D2Fe",], by = "margin")
  #nat_dist and build_ht N.S.

adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+canopy50, data=env[env$tree_ID!="D2Fe",], by = "margin")
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+canopy126, data=env[env$tree_ID!="D2Fe",], by = "margin")

adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+roads50, data=env[env$tree_ID!="D2Fe",], by = "margin") #N.S.
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+roads126, data=env[env$tree_ID!="D2Fe",], by = "margin") #marginally significant

adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+water50, data=env[env$tree_ID!="D2Fe",], by = "margin")
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+water126, data=env[env$tree_ID!="D2Fe",], by = "margin")

adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+built50, data=env[env$tree_ID!="D2Fe",], by = "margin") #marginally significant
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+built126, data=env[env$tree_ID!="D2Fe",], by = "margin") #N.S.

adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+roads126+built50+canopy50, data=env[env$tree_ID!="D2Fe",], by = "margin") #N.S.
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+roads126+built50+canopy126, data=env[env$tree_ID!="D2Fe",], by = "margin")

adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+roads126+built50+canopy126+water50, data=env[env$tree_ID!="D2Fe",], by = "margin")
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+roads126+built50+canopy126+water126, data=env[env$tree_ID!="D2Fe",], by = "margin")
  #roads126 and water126 marginally significant
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+built50+canopy126+water50, data=env[env$tree_ID!="D2Fe",], by = "margin")
  #built50 marginally significant
adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+canopy126+water50, data=env[env$tree_ID!="D2Fe",], by = "margin")

write.csv(
  rbind(
    adonis(nfvfnf.bray$bray~mistletoe+tree_ht+build_ht+canopy126+tree_type, data=env[env$tree_ID!="D2Fe",], by = "margin")$aov.tab,
    adonis(nfvfnf.bray$bray.bal~mistletoe+tree_ht+build_ht+canopy126+tree_type, data=env[env$tree_ID!="D2Fe",], by = "margin")$aov.tab,
    adonis(nfvfnf.bray$bray.gra~mistletoe+tree_ht+build_ht+canopy126+tree_type, data=env[env$tree_ID!="D2Fe",], by = "margin")$aov.tab,
    adonis(nfvfnf.beta$beta.sor~mistletoe+tree_ht+canopy126+water50+tree_type, data=env[env$tree_ID!="D2Fe",], by = "margin")$aov.tab,
    adonis(nfvfnf.beta$beta.sim~mistletoe+tree_ht+canopy126+water50+tree_type, data=env[env$tree_ID!="D2Fe",], by = "margin")$aov.tab,
    adonis(nfvfnf.beta$beta.sne~mistletoe+tree_ht+canopy126+water50+tree_type, data=env[env$tree_ID!="D2Fe",], by = "margin")$aov.tab
  ), "permanova_nfvfnf.csv", quote=FALSE
)

## fignf vs. figf

fnfvff.bray<-bray.part(fnfvff.com)

fnfvff.com2<-ifelse(fnfvff.com>0,	1, 0)
fnfvff.beta<-beta.pair(fnfvff.com2)

only1.tab<-colSums(table(figging, tree_ID2))
tree_ID3<-replace(tree_ID2, tree_ID2 %in% names(only1.tab)[only1.tab==1], "none")

write.csv(
  rbind(
    #with abundances
    
    adonis(fnfvff.bray$bray~figging+tree_sp2, permutations=how(within = Within(type = "free"), plots = Plots(strata = tree_ID3), nperm=999), by = "margin")$aov,
    adonis(fnfvff.bray$bray.bal~figging+tree_sp2, permutations=how(within = Within(type = "free"), plots = Plots(strata = tree_ID3), nperm=999), by = "margin")$aov,
    adonis(fnfvff.bray$bray.gra~figging+tree_sp2, permutations=how(within = Within(type = "free"), plots = Plots(strata = tree_ID3), nperm=999), by = "margin")$aov,
    
    #presence-absence only
    
    adonis(fnfvff.beta$beta.sor~figging+tree_sp2, permutations=how(within = Within(type = "free"), plots = Plots(strata = tree_ID3), nperm=999), by = "margin")$aov,
    adonis(fnfvff.beta$beta.sim~figging+tree_sp2, permutations=how(within = Within(type = "free"), plots = Plots(strata = tree_ID3), nperm=999), by = "margin")$aov,
    adonis(fnfvff.beta$beta.sne~figging+tree_sp2, permutations=how(within = Within(type = "free"), plots = Plots(strata = tree_ID3), nperm=999), by = "margin")$aov
    
  ), "permanova_fnfvff.csv", quote=FALSE
)

## Beta dispersion

all.dist<-vegdist(all.com, method = "bray")

all.betadisper<-betadisper(all.dist, tree_type2)
permutest(all.betadisper)
TukeyHSD(all.betadisper)

all.dist2<-vegdist(all.com2, method = "bray")

all.betadisper2<-betadisper(all.dist2, tree_type2)
permutest(all.betadisper2)
TukeyHSD(all.betadisper2)

## figures

col.tree_type<-ifelse(env$tree_type=="nonfig", "blue", "darkgreen")
pch.tree_type<-ifelse(env$tree_type=="nonfig", 17, 16)

col.figging<-ifelse(figging=="figf","red","darkgreen")
pch.figging<-ifelse(figging=="figf",1,16)

pch.tree_type2<-c(pch.tree_type[env$tree_ID!="D2Fe"], pch.figging[figging=="figf"])
col.tree_type2<-c(col.tree_type[env$tree_ID!="D2Fe"], col.figging[figging=="figf"])

tree_type2<-c(as.character(with(env, tree_type[tree_ID!="D2Fe"])), as.character(figging[figging=="figf"]))
tree_type2<-replace(tree_type2, tree_type2=="fig", "fignf")

### with abundances

all.com<-rbind(nfvfnf.com[env$tree_ID!="D2Fe",], fnfvff.com[figging=="figf",])

all.mds<-metaMDS(all.com, k=3, distance = "bray")

### presence-absences only

all.com2<-ifelse(all.com>0, 1, 0)

all.mds2<-metaMDS(all.com2, k=3, distance = "bray")

windows()
par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.5,1,0), cex.lab=1.5)
plot(all.mds, display="species", type="n", choices=1:2)
points(all.mds, display = "sites", pch=pch.tree_type2, col=col.tree_type2, choices=1:2, cex=2, lwd=2)
text(all.mds, display = "species", col=ifelse(rownames(all.mds$species) %in% rownames(with(nfvfnf.indval, indval[pval<0.05,])) |
                                                rownames(all.mds$species) %in% rownames(with(fnfvff.indval, indval[pval<0.05,])),
                                              "black", "gray"), choices=1:2, cex=0.7)
legend("topleft", legend = "(a)", inset = c(-0.4, -0.1), xpd = TRUE, bty = "n", cex=1.5)
plot(all.mds, display="species", type="n", choices=c(3,2))
points(all.mds, display = "sites", pch=pch.tree_type2, col=col.tree_type2, choices=c(3,2), cex=2, lwd=2)
text(all.mds, display = "species", col=ifelse(rownames(all.mds$species) %in% rownames(with(nfvfnf.indval, indval[pval<0.05,])) |
                                                rownames(all.mds$species) %in% rownames(with(fnfvff.indval, indval[pval<0.05,])),
                                              "black", "gray"), choices=c(3,2), cex=0.7)
legend("topleft", legend = "(b)", inset = c(-0.4, -0.1), xpd = TRUE, bty = "n", cex=1.5)
plot(all.mds2, display="species", type="n", choices=1:2)
points(all.mds2, display = "sites", pch=pch.tree_type2, col=col.tree_type2, choices=1:2, cex=2, lwd=2)
text(all.mds2, display = "species", col=ifelse(rownames(all.mds$species) %in% rownames(with(nfvfnf.indval, indval[pval<0.05,])) |
                                                 rownames(all.mds$species) %in% rownames(with(fnfvff.indval, indval[pval<0.05,])),
                                               "black", "gray"), choices=1:2, cex=0.7)
legend("topleft", legend = "(c)", inset = c(-0.4, -0.1), xpd = TRUE, bty = "n", cex=1.5)
plot(all.mds2, display="species", type="n", choices=c(3,2))
points(all.mds2, display = "sites", pch=pch.tree_type2, col=col.tree_type2, choices=c(3,2), cex=2, lwd=2)
text(all.mds2, display = "species", col=ifelse(rownames(all.mds$species) %in% rownames(with(nfvfnf.indval, indval[pval<0.05,])) |
                                                 rownames(all.mds$species) %in% rownames(with(fnfvff.indval, indval[pval<0.05,])),
                                               "black", "gray"), choices=c(3,2), cex=0.7)
legend("topleft", legend = "(d)", inset = c(-0.4, -0.1), xpd = TRUE, bty = "n", cex=1.5)

# testing spatial autocorrelation

dist<-as.matrix(dist(coordinates(trees.ly)[env$tree_ID!="D2Fe",1:2], diag=TRUE, upper=TRUE))
dist1<-1/dist
diag(dist1)<-0

c(Moran.I(residuals(nfvfnf.D0.mod.top1[[1]]), dist1)$p.value,
  Moran.I(residuals(nfvfnf.D1.mod.top1[[1]]), dist1)$p.value,
  Moran.I(residuals(nfvfnf.D2.mod.top1[[1]]), dist1)$p.value,
  Moran.I(residuals(nfvfnf.eabund.mod.top1[[1]]), dist1)$p.value,
  Moran.I(residuals(nfvfnf.Iabund.mod.top1[[1]]), dist1)$p.value)

mantel(nfvfnf.bray$bray, vegdist(coordinates(trees.ly)[env$tree_ID!="D2Fe",1:2], method = "euclidean"))
mantel(nfvfnf.beta$beta.sor, vegdist(coordinates(trees.ly)[env$tree_ID!="D2Fe",1:2], method = "euclidean"))

# iNExt rarefaction

bird.com.sum<-apply(all.com, 2, function(x)
  tapply(x, tree_type2, sum))
bird.com.inext<-apply(bird.com.sum, 1, function(x) sort(x[x!=0], decreasing = TRUE))

bird.inext0<-iNEXT(bird.com.inext, q = 0, endpoint=max(rowSums(bird.com.sum)))

df<-fortify.iNEXT(bird.inext0, type=1)
df.point <- df[which(df$method=="observed"),]
df.line <- df[which(df$method!="extrapolated"),]

windows()
ggplot(df.line, aes(x=x, y=y)) + 
  geom_point(aes(shape=site, colour=site), size=5, stroke=3, data=df.point) +
  geom_line(aes(linetype=method, colour=site), lwd=1.5, data=df.line, show.legend = FALSE) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=site, colour=NULL), alpha=0.2, show.legend = FALSE) +
    theme_classic()+
  theme(legend.position = c(0.8,0.2), legend.key.size = unit(1, units="cm"),
        text=element_text(size=18),
        axis.title.y = element_text(margin=margin(0,10,0,0)))+
  labs(x="Number of individuals", y="Species richness", shape="Tree Type", colour="Tree Type") + 
  scale_shape_manual(values=c(1,16,17), labels=c("FIGF", "FIGNF", "NONFIG"))+
  scale_colour_discrete(labels=c("FIGF", "FIGNF", "NONFIG"))

# indicator species

nfvfnf.indval<-indval(nfvfnf.com[env$tree_ID!="D2Fe",colSums(nfvfnf.com)!=0], env$tree_type[env$tree_ID!="D2Fe"])

with(nfvfnf.indval, indval[pval<0.05,])
with(nfvfnf.indval, relfrq[pval<0.05,])
with(nfvfnf.indval, relabu[pval<0.05,])

fnfvff.indval<-indval(fnfvff.com[,colSums(fnfvff.com)!=0], figging)

with(fnfvff.indval, indval[pval<0.05,])
with(fnfvff.indval, relfrq[pval<0.05,])
with(fnfvff.indval, relabu[pval<0.05,])

mistletoe.indval<-indval(nfvfnf.com[env$tree_ID!="D2Fe",colSums(nfvfnf.com)!=0], env$mistletoe[env$tree_ID!="D2Fe"])

with(mistletoe.indval, indval[pval<0.05,])
with(mistletoe.indval, relfrq[pval<0.05,])
with(mistletoe.indval, relabu[pval<0.05,])

nfvfnf.I.indval<-indval(nfvfnf.I.com[env$tree_ID!="D2Fe",colSums(nfvfnf.I.com)!=0], env$tree_type[env$tree_ID!="D2Fe"])

with(nfvfnf.I.indval, indval[pval<0.05,])
with(nfvfnf.I.indval, relfrq[pval<0.05,])
with(nfvfnf.I.indval, relabu[pval<0.05,])

fnfvff.I.indval<-indval(fnfvff.I.com[,colSums(fnfvff.I.com)!=0], figging)

with(fnfvff.I.indval, indval[pval<0.05,])
with(fnfvff.I.indval, relfrq[pval<0.05,])
with(fnfvff.I.indval, relabu[pval<0.05,])