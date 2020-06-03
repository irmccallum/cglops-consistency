#-----------------------------
# validation_brick.R
# validate CGLOPS 100m Africa map
# with Copernicus HotSpot maps
# IIASA, Apr 28 2020
#-----------------------------

library(raster)
library(plyr)
library(caret)
library(gdalUtils)
library(ggpubr)
library(reshape)

st.out <- NULL
cf.out <- NULL
m.list <- NULL

# load merged hotspots and legend
s <- shapefile("LC_merged.shp")
s.l <- read.csv("hotspots_rev.csv")
s <- merge(s, s.l, by='class_name')

# load LC and legend
r <- raster("ProbaV_LC100_epoch2015-base_Africa_v2.1.1_discrete-classification_EPSG-4326.tif")
names(r) <- "code"
r.l <- read.csv("lccs_new_id.csv")

for (i in 1: length(unique(s$aoi_name))){

# select one hotspot to work with
index <- s$aoi_name == unique(s$aoi_name)[i] 
s.map <- s[index,]
shapefile(s.map, "tempshp.shp", overwrite = T)

print(unique(s$aoi_name)[i])

# crop raster
c <- crop(r,s.map)
c <- mask(c,s.map)

# add legend to raster
c <- subs(c, r.l, by='code', which='LC_ID')

# rasterize hotspot maps
gdal_rasterize("tempshp.shp","temphs.tif", tap = TRUE,a="hotspid",l="tempshp", a_nodata=0, verbose = TRUE,te=c(xmin(s.map),ymin(s.map),xmax(s.map),ymax(s.map)),tr=c(res(r)[1], res(r)[2]))
sr <- raster("temphs.tif")
writeRaster(c,"ctemp.tif", overwrite = TRUE )
align_rasters("temphs.tif","ctemp.tif","align.tif", overwrite = TRUE)
sr <- raster("align.tif")

# stack rasters
st <- brick(c, sr)
names(st) <- c("LC100","Hotspots")

lc.df = as.data.frame(st,xy=TRUE)

mdata <- melt(lc.df, id=c("x","y"))
mdata$value <- as.factor(mdata$value)
levels(mdata$value) <- c(1:7)

# plot land cover maps
p1 <- ggplot() + 
  geom_raster(data=mdata,aes(x,y,fill=value)) +
  facet_wrap(~ variable)+
  scale_fill_manual(values = c('orange', 'green', 'pink',"red","blue","brown","yellow"),
                    labels = c('Cultivated & Managed', 'Natural Vegetation', 'Wetlands',"Urban","Water","Bare","Snow & Ice"),
                    breaks = c("1","2","3","4","5","6","7"),
                    limits = c("1","2","3","4","5","6","7")) +
  theme(aspect.ratio = 1) +
  ggtitle(unique(s$aoi_name)[i])

ggsave(paste0("./results/",unique(s$aoi_name)[i],"LC.png",sep=""),  width=10, height=6, dpi=300)

# agreement maps (agreement = 0; disagreement = 1)
agm <- st$LC100 - st$Hotspots
agm[agm != 0] <- 1
agm.df = as.data.frame(agm,xy=TRUE)

# plot agreement maps
p <- ggplot()+ 
  geom_raster(data=agm.df,aes(x,y,fill=layer)) +
  theme(aspect.ratio = 1) +
  scale_fill_viridis_c(option = "cividis",na.value = "transparent",breaks = c(0,1),labels = c("agree", "disagree")) +
  ggtitle(unique(s$aoi_name)[i]) 
ggsave(paste0("./results/",unique(s$aoi_name)[i],"AGREE.png",sep=""),  width=5, height=5, dpi=300)

# convert stack to dataframe
st.df <- as.data.frame(st)
names(st.df) <- c("ID_LC","ID_hs")
st.df <- na.omit(st.df)

l.hs <- read.csv("legend_hotspots.csv")
l.lc <- read.csv("legend_lc100.csv")

# join class names to df
st.df <- join(st.df,l.lc, by = "ID_LC")
st.df <- join(st.df,l.hs, by = "ID_hs")
st.df <- st.df[,3:4]

# write out stats
cf <- caret::confusionMatrix(st.df$LC100, st.df$Hotspots)
cf.out <- rbind(cf.out,cf$overall)
write.csv(cf.out,file=paste0("./results/ALLSTATS.csv"))

# write out matrix
st.tab <- table(st.df$LC100,st.df$Hotspots)
st.out <- rbind(st.out,st.tab)
write.csv(st.out,paste0("./results/ALLTABLES.csv"))

}
