# SETUP ####

# load packages
library(tidyverse)
library(phyloseq)
library(patchwork)
library(vegan)
library(gdm)


# functions
source("./R/functions.R")

# options
theme_set(theme_bw() +
            theme(strip.text = element_text(size=18,face='bold'),
                  axis.text.x = element_text(face='bold',size=14,angle=90,hjust=1,vjust=.5),
                  axis.text.y = element_text(face='bold',size=14),
                  axis.title = element_text(face='bold',size=18),
                  strip.background = element_blank())
)

# physeq data
fungi_spp <- readRDS("./output/fungal_species_full_phyloseq_object.RDS")


ord <- ordinate(fungi_spp,method = "NMDS")
data.frame(MDS1=ord$points[,1],
           MDS2=ord$points[,2],
           size_fraction = fungi_spp@sam_data$size_fraction_um) %>% 
  ggplot(aes(x=MDS1,y=MDS2,color=size_fraction)) +
  geom_point(size=3) +
  scale_color_viridis_d(option='turbo',end=.8,begin=.2) +
  labs(color="Size fraction") +
  theme(legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12))
ggsave("./output/figs/nmds_plot_size_fraction.png",dpi=300,height = 6,width = 6)

# add bullshit long and lat
fungi_spp@sam_data$latitude <- 22.0854
fungi_spp@sam_data$longitude <- -157.9821
fungi_melt <- psmelt(fungi_spp)



# GDM ####
# extract species by site info
fungi_ra_melt <- 
  fungi_spp %>% 
  transform_sample_counts(function(x){ifelse(x>0,1,0)}) %>% 
  psmelt()


# biological data
# get columns with xy, site ID, and species data
sppTab <- fungi_ra_melt %>% dplyr::select(OTU,Sample,size_fraction_um,depth_m,Abundance)

# get columns with site ID, env. data, and xy-coordinates
envTab <- fungi_ra_melt %>% dplyr::select(Sample, longitude, latitude,
                                          all_of(c("ug_n","d15n_vs_air","ug_c","d13c_pc_vs_vpdb","ug_poc",
                                                          "d13c_poc_vs_vpdb","l_filtered","ug_n_l","ug_c_l","ug_poc_l")))

# remove sites with missing envir data
sppTab <- 
  sppTab[complete.cases(envTab),]
envTab <- 
  envTab[complete.cases(envTab),]


# format for gdm
gdmTab <- formatsitepair(bioData=sppTab, 
                         bioFormat=2, #x-y spp list
                         XColumn="longitude", 
                         YColumn="latitude",
                         sppColumn="OTU", 
                         siteColumn="Sample", 
                         predData=envTab,
                         abundance = TRUE,
                         abundColumn = "Abundance")

# fit GDM
gdm <- gdm(data = gdmTab,geo = FALSE)

# quick look at model fits
summary(gdm)

# predictions from model (using same distances)
gdm_pred <- predict(object=gdm, data=gdmTab)
preds <- data.frame(observed = gdmTab$distance,
                    predicted = gdm_pred,
                    dist = gdm$ecological,gdmTab)
p1 <- preds %>% 
  ggplot(aes(x=observed,y=predicted)) +
  geom_point() +
  geom_smooth() + ggtitle("Overall")

# Extract splines
isplines <- gdm::isplineExtract(gdm) %>% as.data.frame()
# rename
names(isplines) <- 
  names(isplines) %>% 
  str_replace("x.","actual_") %>% 
  str_replace("y.","partial_")

# important variables: ug_n | d13c_pc_vs_vpdb | ug_poc 
# plot
p1 <- 
  isplines %>% 
  ggplot(aes(x=actual_ug_n,y=partial_ug_n)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75) +
  labs(y="f(ug N)",x="ug N")

p2 <- 
  isplines %>% 
  ggplot(aes(x=actual_d13c_pc_vs_vpdb,y=partial_d13c_pc_vs_vpdb)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75) +
  labs(y="f(d13c_pc_vs_vpdb)",x="d13c_pc_vs_vpdb")

p3 <- 
  isplines %>% 
  ggplot(aes(x=actual_ug_poc,y=partial_ug_poc)) +
  # geom_vline(xintercept = 100,linetype=2,color='gray') +
  geom_smooth(se=FALSE,linewidth=2,alpha=.75) +
  labs(y="f(ug POC)",x="ug POC")
p1 | p2 | p3
ggsave("./output/figs/GDM_Splines.png",width = 12,height = 8)
  