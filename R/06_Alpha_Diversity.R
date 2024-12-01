# SETUP ####

# load packages
library(tidyverse)
library(phyloseq)
library(patchwork)
library(vegan)

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
fungi_melt <- psmelt(fungi_spp)


# fungi_spp@sam_data %>% View
# fungi_spp@otu_table %>% View


# ALPHA DIV ####

# build data frame
alpha_df <- 
  fungi_spp@sam_data %>% as("data.frame") %>% 
  mutate(richness = estimate_richness(fungi_spp,measures = "Observed") %>% 
           pluck("Observed"))

# plot
alpha_df %>% 
  ggplot(aes(x=depth_m,y=richness,color=depth_m)) +
  geom_point(size=3) +
  facet_wrap(~size_fraction_um,scales = 'free')

# regression
alpha_df %>% 
  glm(data=.,
      formula=richness ~ depth_m + size_fraction_um) %>% 
  broom::tidy() %>% 
  saveRDS("output/alpha_div_glm_table.RDS")
  # no sig diff in richness based on size fraction or depth


# CONNECTIVITY ####
# community membership overlap across depths

s <- 
fungi_spp %>% 
  subset_samples(size_fraction_um == "0.3-1")
s <- 
  s %>% 
  subset_taxa(taxa_sums(s) > 0)
m <- 
  fungi_spp %>% 
  subset_samples(size_fraction_um == "1-53")
m <- 
  m %>% 
  subset_taxa(taxa_sums(s) > 0)
l <- 
  fungi_spp %>% 
  subset_samples(size_fraction_um == ">53")
l <- 
  l %>% 
  subset_taxa(taxa_sums(s) > 0)

# store boolean for taxa present in all size fractions (in kingdom slot)
fungi_spp@tax_table[,1] <- 
taxa_names(fungi_spp) %in% taxa_names(s) &
  taxa_names(fungi_spp) %in% taxa_names(m) &
  taxa_names(fungi_spp) %in% taxa_names(l)

merg <- fungi_spp %>%
  merge_samples("size_fraction_um")
merg@sam_data$size_fraction_um <- sample_names(merg) %>% factor(levels=c("0.3-1", "1-53", ">53"))
merg %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Class",x="size_fraction_um") +
  facet_wrap(~Kingdom) +
  scale_fill_viridis_d(option = "turbo") +
  ggtitle("Ubiquity across filter sizes") +
  labs(caption = "The smallest filter size captured 100% of the detected fungi.\nLarger filter sizes failed to capture some fungi.",
       x="Size fraction (um)",y="Relative abundance") +
  theme(legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12))
ggsave("./output/figs/fungal_ubiquity_across_filter_size.png",dpi=300,height = 8,width = 8)
