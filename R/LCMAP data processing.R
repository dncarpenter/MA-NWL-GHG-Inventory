## Load packages
library(terra)
library(tidyverse)
# library(tigris)
library(sf)

## Set local data directories for raw LCMAP data and processed/output (dif. from vector data which doesn't take up much space)
env.r <- "C:/Users/Dunbar.Carpenter/OneDrive - Commonwealth of Massachusetts/Data & Tools/LULC/LCMAP/Data/"
env.p <- "C:/Users/Dunbar.Carpenter/OneDrive - Commonwealth of Massachusetts/Analyses/Mass NWL GHG Inventory/data/"
aws.r <- "G:/LCMAP_Data/raw/"
aws.p <- "G:/LCMAP_Data/"

dat.raw.dir = env.r
dat.prc.dir = env.p


#### LCMAP
## Load 2016 primary & secondary land cover data
lcmap1 <- rast(paste(dat.raw.dir, 'LCMAP_CU_2016_V12_LCPRI.tiff', sep=''))
lcmap2 <- rast(paste(dat.raw.dir, 'LCMAP_CU_2016_V12_LCSEC.tiff', sep=''))
lcmap1 <- rast(paste(dat.raw.dir, 'LCMAP_CU_2020_V12_LCPRI.tiff', sep=''))
lcmap2 <- rast(paste(dat.raw.dir, 'LCMAP_CU_2020_V12_LCSEC.tiff', sep=''))

lcmap <- c(lcmap1, lcmap2)

## Set levels for categorical rasters
lcmap.lgd <- data.frame(id=1:8, cover=c('Developed','Cropland','Grass/Shrub','Tree Cover',
                                          'Water','Wetland','Ice/Snow','Barren'), value=1:8)
levels(lcmap) <- list(lcmap.lgd, lcmap.lgd)
names(lcmap) <- c('cover1','cover2')


## AOI boundary: Load geometry, re-project, convert to spat vector, make 500m buffer
# ma.bdry <- states(cb = T, year = 2020) %>% 
#                   filter_state("Massachusetts") %>%
aoi.bdry <- st_read(dsn = "input data/OUTLINE25K_POLY.shp") %>%
                      st_geometry() %>%
                      st_transform(crs = crs(lcmap)) 
aoi <- vect(aoi.bdry) # make spatvector v. of AOI bdry
aoi.bfr <- st_buffer(aoi.bdry, dist = 500) %>% # add buffer & make spatvector from AOI bdry
                    st_union() %>%
                    vect()

## Crop LCMAP data to buffered AOI boundary
aoi.lc <- mask(crop(lcmap, aoi.bfr), mask = aoi.bfr)

## plot
plot(aoi.lc$cover1, legend='topright')
plot(aoi, add=T)

## Make land cover frequency table
lc.frq <-  freq(aoi.lc$cover1) # 
lc.frqtbl <- tibble(ClassScheme=rep('Original',nrow(lc.frq)),
                      Class=lc.frq$value,
                      Area.ha=round(lc.frq[,'count']*res(aoi.lc)[1]^2/100^2),
                      Area.ac=round(lc.frq[,'count']*res(aoi.lc)[1]^2/4046.86),
                      Area.sqkm=round(lc.frq[,'count']*res(aoi.lc)[1]^2/1000^2),
                      Area.sqmi=round(lc.frq[,'count']*res(aoi.lc)[1]^2/2590000),
                      Area.pct=round(100*lc.frq[,'count']/sum(lc.frq[,'count']),2)) %>%
  arrange(desc(Area.ha))
# lc.frqtbl

#### I. Forested Wetland Modification
## Make LCPRI wetland, LCSEC forest, and combo rasters; plot
malc1.wl <- classify(aoi.lc$cover1, rcl=cbind(6,6), othersNA=T)
malc2.tr <- classify(aoi.lc$cover2, rcl=cbind(4,4), othersNA=T)
lc.forwet <- malc1.wl*malc2.tr

# plot(malc1.wl, col='blue1')
# plot(malc2.tr, add=T, col='green1')
# plot(lc.forwet, add=T, col='red3')

## Add modified LC layer to raster stack
lc.forwet <- subst(lc.forwet, from=NA, to=0)
aoi.lc$cover3 <- aoi.lc$cover1 + lc.forwet
aoi.lc$cover3 <- subst(aoi.lc$cover3, from=30, to=4)
aoi.lc3.lgd <- lcmap.lgd %>%
  rename('cover3'='cover')
levels(aoi.lc$cover3) <- aoi.lc3.lgd ## apply legend from original LC data
coltab(aoi.lc$cover3) <- coltab(aoi.lc$cover1) ## apply color table from original LC data
# plot(aoi.lc)
# plot(aoi.lc$cover3, legend='topright')
writeRaster(aoi.lc, 'C:/Users/Dunbar.Carpenter/OneDrive - Commonwealth of Massachusetts/Analyses/Mass NWL GHG Inventory Processing/data/intermediate/LMCAP MA landcover 2020 wetland-forest adjusted.tif', overwrite=T)
# LMCAP MA landcover 2016 wetland-forest adjusted.tif

## Append frequency table for new LC layer to original for comparison
malc3.frq <- freq(mask(aoi.lc$cover3, ma))
lc.frqtbl <- tibble(ClassScheme=rep('ForWet mod',nrow(malc3.frq)),
                      Class=malc3.frq$value,
                      Area.ha=round(malc3.frq[,'count']*res(aoi.lc)[1]^2/100^2),
                      Area.ac=round(malc3.frq[,'count']*res(aoi.lc)[1]^2/4046.86),
                      Area.sqkm=round(malc3.frq[,'count']*res(aoi.lc)[1]^2/1000^2),
                      Area.sqmi=round(malc3.frq[,'count']*res(aoi.lc)[1]^2/2590000),
                      Area.pct=round(100*malc3.frq[,'count']/sum(malc3.frq[,'count']),2)) %>%
  arrange(desc(Area.ha)) %>%
  bind_rows(lc.frqtbl) #%>%
# write_csv('Analyses/Mass LULC/data/intermediate/LCMAP2016 MA LC scheme comparison.csv')
# lc.frqtbl


############ II. Settlement Land
#### small set extent & crop MA LC data
plot(aoi.lc$cover3, legend = 'topright')
# sa.e <- draw()  ## manually select sample area
# e <- ext(aoi.lc)
e <- ext(2014000, 2067000, 2354000, 2397000) ## larger sample area ext (~10% of state)
# e <- ext(1977495, 1981335, 2394585, 2397105) ## smaller sample area ext (~0.04% of state)

lc <- crop(aoi.lc$cover3, ext(e)) ## crop to relevant extent

# (aoi.lc.area <- sum(freq(aoi.lc$cover3)$count)*30^2/1000^2) ## km^2 area of LC data for all MA
# (lc.area <- sum(freq(lc)$count)*30^2/1000^2)              ## km^2 area of LC data  for extent e
# 100*lc.area/aoi.lc.area                                    ## sample area as % of MA
plot(lc)


## A. LCPRI = Developed
aoi.lc$cover4 <- aoi.lc$cover3
levels(aoi.lc$cover4) <- data.frame(id=1:8, class=c('Settlement','Cropland','Grass/Shrub','Tree Cover',
                                                   'Water','Wetland','Ice/Snow','Barren'), value=1:8)

##	B. NWL Patches in Settlement Land
#####
## IF LCPRI in {Cropland, Grass/Shrub, Wetland, Water, Barren, Tree Cover} 
## AND patch area <= 4 ha (~7x7 30m pixels) 
## AND surrounded by LCPRI Developed

## make no-settlement mask, assigning 1 to all non-settlement cells
ns <- classify(lc, cbind(0:8, c(NA, NA,rep(1,7))))
plot(ns)

system.time({
  ns.ptchs <- patches(ns, directions=4)
})
# plot(ns.ptchs, main='Patch ID')
## freqency table 
ns.ptch.frq <- freq(ns.ptchs) %>%
  as_tibble()
# ns.ptch.frq

## non-settlement patch size
ns.psz <- classify(ns.ptchs, ns.ptch.frq[,-1])

## calculate sieve size (~4ha)
ns.siv.sz <- round(sqrt(4*100^2/prod(res(lc))))^2 ## no of patches <= sieve size # of cells
# length(which(ns.ptch.frq$count <= ns.siv.sz))

## values (cell nos) for patches <= X cells in area (X=sieve size)
ns.sm.ptch <- ns.ptch.frq$value[which(ns.ptch.frq$count <= ns.siv.sz)]
## New non-settlement mask to be sieved; assign new value to all cells in patches <= sieve size
ns.siv <- ns
ns.siv[ns.ptchs %in% ns.sm.ptch] <- 2
#######
#### Examine sieve results
########
# freq(ns.siv)
plot(lc)
plot(ns.siv, main='Non-settlement Patch Sieve')

## Select area to zoom in on to a small extent to check the results
# e <- draw() 
e2 <- ext(lc)
plot(ns.ptchs, main='Patch ID', ext=e2)
plot(ns.psz, col=hcl.colors(nrow(ns.ptch.frq), palette='Spectral'), main='Patch Size', ext=e2)
plot(ns.siv, main='Non-settlement Patch Sieve', ext=e2)
plot(lc, ext=e2)


########
## C. Tree cover
## IF LCPRI == Tree Cover
## AND patch area <= 100 ha (~7x7 30m pixels) 
## AND surrounded by LCPRI Developed
#######
tc <- segregate(lc, classes=4, other=NA)

## ID tree cover patches
system.time({
  tc.ptchs <- patches(tc, directions=4, allowGaps=F)
})
# plot(tc.ptchs, main='Patch ID')
## freqency table 
tc.ptch.frq <- freq(tc.ptchs) %>%
  as_tibble()
# tc.ptch.frq

## tree cover patch size
tc.psz <- classify(tc.ptchs, tc.ptch.frq[,-1])
plot(tc.psz, col=hcl.colors(nrow(tc.ptch.frq), palette='Spectral'), main='Patch Size', ext=e2)

## sieve area (ha), used to calculate sieve size in cells 
tc.siv.ha = 40
tc.siv.cells <- round(sqrt(tc.siv.ha*100^2/prod(res(lc))))^2 
# length(which(tc.ptch.frq$count <= tc.siv.ha)) ## no of patches <= sieve size # of cells

## values (cell nos) for small patches <= sieve size (in cells)
tc.smp.vals <- tc.ptch.frq$value[which(tc.ptch.frq$count <= tc.siv.cells)]
## New non-settlement mask to be sieved; assign new value to all cells in patches <= sieve size
tc.siv <- tc
tc.siv[tc.ptchs %in% tc.smp.vals] <- 2
# freq(tc.siv)


plot(lc, main='Land Cover')
plot(tc, main='Tree Cover', type='classes', levels=c('tree cover'))
plot(tc.siv, main='Big & Small Tree Cover Patches', type='classes', 
     levels=c(paste('Tree cover, \npatches >', tc.siv.ha, 'ha'), paste('Tree cover, \npatches <=', tc.siv.ha, 'ha')))
# plot(tc.sm.ptchs, main='Small Tree Cover Patches', type='classes', 
#      levels=c(paste('Tree cover, \npatches <=', tc.siv.ha, 'ha')))

## Isolate small forest patches
tc.sp <- subst(tc.siv, from=1, to=NA) # small tree cover patches only
stl.nstl <- classify(lc, cbind(c(NA, 1:8), c(0, 1, rep(2,7)))) # classify lc as settlement v. non-settlement land
sn.stp.msk <- mask(stl.nstl, mask = tc.sp, maskvalues=2) # mask small tree patches from settlement v. non-settlement raster

plot(lc, main='Land Cover')
plot(stl.nstl, main='Settlement v. Non-Settlement', type='classes', levels=c('settlement','non-settlement'))
plot(tc.sp, main='Small Tree Cover Patches', type='classes',
     levels=c( paste('Tree cover, \npatches <=', tc.siv.ha, 'ha')))
plot(sn.stp.msk, main=paste('Settlement + Non-Settlement - Small Tree Patches (<=', tc.siv.ha, 'ha)'), 
     type='classes', levels=c('settlement','non-settlement'))

##  Focal modal reclassification: calculate most frequent non-NA value w/in 3x3 window of all NA cells 
## (where NA=small tree cover patches, 1=settlement, 2=non-settlement, 0=water)
sn.stpm.fm <- focal(sn.stp.msk, w=3, fun='modal', na.rm=T, na.policy='only')
## show only focal modal results i.e. most frequent neighbor of cells in small tree cover patches;
stp.fm <- mask(sn.stpm.fm, mask = tc.sp)

plot(sn.stpm.fm, main='Settlement/Non-Settlement + Small Tree Patches Edges Reclassified', 
     type='classes', levels=c('settlement','non-settlement'))
plot(stp.fm, main='Reclassified Small Tree Patch Edges', type='classes', 
     levels=c('settlement \nreclass','non-settlement \nreclass'))
plot(lc, main='Land Cover')

##
stp.bdry <- boundaries(sn.stp.msk, inner=F, directions=8, falseval=NA) # get small tree patch boundaries
stp.pid <- mask(tc.ptchs, mask = tc.sp) 
stp.pid.bdry <- mask(tc.ptchs, mask = stp.bdry) # patch IDs for all small tree patch boundary cells


plot(stp.bdry, main='Small Tree Patch Boundaries', type='classes', levels=c('Small tree\npatch boundary'))
plot(stp.pid, main='Small Tree Patch IDs')
plot(stp.pid.bdry, main='Patch IDs for Small Tree Patch Boundaries')


## get modal neighbor class cell counts for each small tree cover patch w/ crosstab; rearrange to wide format
stp.fm.rc <- subst(stp.fm, from=NA, to=0) # convert NA to 0
stp.xtb <- crosstab(c(stp.pid.bdry, stp.fm), long=T, useNA=T)  %>% 
  as_tibble() %>%
  pivot_wider(names_from = focal_modal, values_from = Freq, values_fill = 0) %>%
  drop_na(patches)
## add column with the most frequent modal neighbor class in each patch
stp.xtb %<>%
  select(`1`,`2`,`NA`) %>%                # select only cell class cols
  apply(MARGIN = 1, FUN = function(x) { if(x[1] > sum(x)/2){return(1)} else {return(2)} }) %>%  # which.max# which col has max cell count? ties return 1st true val, so 2 (non-settlement) or NA
  # c(2, NA, 1)[.] %>%                      # get names of these cols (`2`=non-settlement, `NA`=NA, `1`=settlement)
  mutate(stp.xtb, MaxNbrCls = .)          # add new col 'MaxNbrCls' with majority surrounding class for each patch
## reclass small tree patches by focal modal 
stp.rc <- classify(stp.pid, rcl = stp.xtb[,c('patches','MaxNbrCls')])  

plot(stp.rc, main='Reclassified Small Tree Patches', type='classes', 
     levels=c('Reclassified \nas settlement','Remains \ntree cover'))
plot(lc, main='Land Cover')
# text(stp.rc, cex=.5)

#######

## Add modified LC layer to raster stack
stl <- classify(lc, cbind(1:8, c(NA, rep(1,7))))
stl2 <-
  malc.forwet <- subst(malc.forwet, from=NA, to=0)
aoi.lc$cover3 <- aoi.lc$cover1 + malc.forwet
aoi.lc$cover3 <- subst(aoi.lc$cover3, from=30, to=4)
aoi.lc3.lgd <- lcmap1.lgd %>%
  rename('cover3'='cover1')
levels(aoi.lc$cover3) <- aoi.lc3.lgd ## apply legend from original LC data
coltab(aoi.lc$cover3) <- coltab(aoi.lc$cover1) ## apply color table from original LC data
# plot(aoi.lc)

########
##	D. All land within 30m (1 pixel) of LCPRI Developed
########
stl <- segregate(lc, classes=1, other=NA)
plot(stl)

stl.bdry <- boundaries(stl, inner=F, directions=4, falseval=NA)

plot(stl.bdry)



############ sandbox
####### Alt. patch analysis w/ landscapemetrics
library(landscapemetrics)
show_patches(sa.ns, directions = 4)
show_lsm(sa.ns, what='lsm_p_area', directions = 4)

system.time({
  ns.pa <- lsm_p_area(sa.ns, directions = 4 )
  ns.ptchs <- get_patches(sa.ns, directions = 4, return_raster = F)
})
ns.ptchs$layer_1$class_1@data@values
########

r1 <- r2 <- rast(ncols=36, nrows=18)
values(r1) <- 1:ncell(r1)
values(r2) <- runif(ncell(r2))
r2 <- classify(r2, cbind(-Inf, 0.5, NA))
r3 <- cover(r2, r1)


p <- vect(system.file("ex/lux.shp", package="terra"))
e <- as.polygons(ext(6, 6.4, 49.75, 50))
values(e) <- data.frame(y=10)
cv <- cover(p, e)
plot(cv, col=rainbow(12))
ci <- cover(p, e, identity=TRUE)
lines(e, lwd=3)

plot(ci, col=rainbow(12))
lines(e, lwd=3)


# create an empty SpatRaster
r <- rast(ncol=10, nrow=10)
# assign values to cells
values(r) <- 1:ncell(r)
s <- r + 10
s <- sqrt(s)
s <- s * r + 5
values(r) <- runif(ncell(r))
r <- round(r)
r <- r == 1

s[r] <- -0.5
s[!r] <- 5
s[s == 5] <- 15

r <- rast(ncol=36,nrow=18)
values(r) <- NA
r[500] <- 1
b <- buffer(r, width=5000000)
plot(r)
plot(b)

set.seed(0)
r <- rast(nrows=10, ncols=10)
values(r) <- sample(3, ncell(r), replace=TRUE)
is.factor(r)

cls <- c("forest", "water", "urban")
# make the raster start at zero
x <- r - 1
levels(x) <- cls
names(x) <- "land cover"
is.factor(x)
x

plot(x, col=c("green", "blue", "light gray"))
text(x, digits=3, cex=.75, halo=TRUE)



lcmap16.frq <- freq(lcmap16.ma)
lcmap.16 <- tibble(Data=rep('LCMS land cover', nrow(lcmap16.frq)),
                   Year=rep(2016,nrow(lcmap16.frq)),
                   ClassScheme=rep('Original',nrow(lcmap16.frq)),
                   Class=c('Settlement','Cropland','Grassland','Forest','Water','Wetland','Other'),
                   Area.ha=round(lcmap16.frq[,'count']*res(lcmap16.ma)[1]^2/100^2),
                   Area.ac=round(lcmap16.frq[,'count']*res(lcmap16.ma)[1]^2/4046.86),
                   Area.sqkm=round(lcmap16.frq[,'count']*res(lcmap16.ma)[1]^2/1000^2),
                   Area.sqmi=round(lcmap16.frq[,'count']*res(lcmap16.ma)[1]^2/2590000),
                   Area.pct=round(100*lcmap16.frq[,'count']/sum(lcmap16.frq[,'count']),2)) %>%
  arrange(desc(Area.ha))
lcmap.16


####
ggplot(lcmap.16, aes(fill = fct_inorder(Class), values = Area.sqmi/10)) +
  geom_waffle(color = "white", n_rows = 20, flip=T) +
  scale_fill_discrete(type = c("#17620E", "#F76B65", "#1E8A88", "#A1FB57", "blue", "#D2D02F", "#969696")) +
  coord_equal() +
  labs(title = "LCMAP Land Cover") +
  theme_void() +
  theme_enhance_waffle()

##########
all.lulc <- bind_rows(lcms.lu, lcmap.lc, nlcd.lc, ccap.lc, malulc)
write.csv(all.lulc, file = 'Analyses/Mass LULC/data/All Land Use Land Cover Summary.csv', row.names = F)
