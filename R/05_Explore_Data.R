# SETUP ####

# load packages
library(tidyverse)
library(phyloseq)
library(patchwork)
library(vegan)

# functions
source("./R/functions.R")

# options
theme_set(theme_minimal() +
            theme(strip.text = element_text(size=18,face='bold'),
                  axis.text = element_text(face='bold',size=14),
                  axis.title = element_text(face='bold',size=18))
  )

# load data, removing obvious garbage (like Kingdom = NA)
ps_its <- readRDS("./output/ITS_physeq_object.RDS") %>% 
  subset_taxa(!is.na(Kingdom)) %>% clean_ps_taxonomy()
ps_ssu <- readRDS("./output/SSU_physeq_object.RDS") %>% 
  subset_taxa(!is.na(Kingdom)) %>% clean_ps_taxonomy()
ps_lsu <- readRDS("./output/LSU_physeq_object.RDS") %>% 
  subset_taxa(!is.na(Kingdom)) %>% clean_ps_taxonomy()


# MERGE DATA ####
  # how much overlap in taxonomic assignments between amplicons?
full <- merge_phyloseq(ps_its,ps_ssu,ps_lsu) %>% 
  subset_taxa(Kingdom != "k__Mitochondrion")

# get sci names
spp <- paste0(full@tax_table[,6]," ",full@tax_table[,7]) %>% str_remove_all(".__") %>% str_replace_all("unclassified","NA")
fam <- full@tax_table[,5] %>% str_remove(".__") %>% unname %>% paste0(" sp.")
ord <- full@tax_table[,4] %>% str_remove(".__") %>% unname %>% paste0(" sp.")
cla <- full@tax_table[,3] %>% str_remove(".__") %>% unname %>% paste0(" sp.")
phy <- full@tax_table[,2] %>% str_remove(".__") %>% unname %>% paste0(" sp.")
kng <- full@tax_table[,1] %>% str_remove(".__") %>% unname %>% paste0(" sp.")

# replace species assignment with new designation
full@tax_table[,7] <- 
ifelse(spp == "NA NA" | spp == "unclassified unclassified",fam,spp) %>% 
  ifelse(. == "NA sp.",ord,.) %>% 
  ifelse(. == "NA sp.",cla,.) %>% 
  ifelse(. == "NA sp.",phy,.) %>% 
  ifelse(. == "NA sp.",kng,.) %>% 
  str_replace(" NA"," sp.")

full <- 
  full %>% clean_ps_taxonomy()

# remove any that are still totally unclassified
# full <- 
# full %>% 
#   subset_taxa(Species != "unclassified sp." & !is.na(full@tax_table[,2]))
full@tax_table[,6] %>% unique %>% unname
full_spp <- 
  full %>% 
  tax_glom("Species",NArm = FALSE)
full_spp@tax_table[,2] %>% unique %>% unname

# collapse to only fungi (glommed at species level)
fungi_spp <- 
  full_spp %>% 
  subset_taxa(Kingdom == "Fungi")
fungi_spp <- 
  fungi_spp %>% 
  subset_taxa(taxa_sums(fungi_spp) > 0)
fungi_spp <- 
  fungi_spp %>% 
  subset_samples(sample_sums(fungi_spp) > 0)
fungi_spp %>% 
  transform_sample_counts(ra) %>% 
  plot_bar(x="size_fraction_um",fill="Class")
saveRDS(fungi_spp,"./output/fungal_species_full_phyloseq_object.RDS")



# Which fungi found by different methods?
ps_its@sam_data$amplicon <- "ITS"
ps_ssu@sam_data$amplicon <- "SSU"
ps_lsu@sam_data$amplicon <- "LSU"
sample_names(ps_its) <- paste0(sample_names(ps_its),"_ITS")
sample_names(ps_ssu) <- paste0(sample_names(ps_ssu),"_SSU")
sample_names(ps_lsu) <- paste0(sample_names(ps_lsu),"_LSU")




full_separate <- merge_phyloseq(ps_its,ps_ssu,ps_lsu)
full_separate@tax_table[,1] %>% unique %>% unname

data.frame(N_ASVs = ntaxa(full),
           N_Samples = nsamples(full),
           N_Fungal_ASVs = full %>% 
             subset_taxa(Kingdom == "Fungi") %>% ntaxa()) %>% 
  saveRDS("./output/ASV_Stats.RDS")

full@otu_table[,1]


full
fungi <- 
full_separate %>% 
  subset_taxa(Kingdom == "Fungi")
fungi <- 
fungi %>% subset_samples(sample_sums(fungi) > 0)

fungi %>% 
  plot_bar(x="size_fraction_um",fill="Phylum")

mat <- fungi@otu_table %>% as('matrix')
mat[mat>0] <- 1
fungi_pa <- fungi
otu <- otu_table(mat,taxa_are_rows = FALSE)
fungi_pa@otu_table <- otu
ord <- 
fungi_pa %>% 
  ordinate(distance = 'jaccard',method = "RDA")

plot_ordination(fungi_pa,ord,color='amplicon')


# find sample arrangements
meta <- full@sam_data %>% as('data.frame')
  #depth
depth_order <- 
  meta %>% 
  arrange(desc(depth_m)) %>% 
  pluck('sample_id')


full %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Kingdom") +
  scale_fill_viridis_d(option='turbo') +
  # scale_x_discrete(limits=depth_order) +
  facet_wrap(~size_fraction_um,scales='free')

# which fungi are found in each size fraction?
full_spp %>% 
  subset_taxa(Kingdom == "Fungi") %>% 
  transform_sample_counts(ra) %>% 
  plot_bar(fill="Class") +
  scale_fill_viridis_d(option='turbo') +
  # scale_x_discrete(limits=depth_order) +
  facet_wrap(~size_fraction_um,scales='free') +
  labs(y="Relative abundance") +
  theme(legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12),
        axis.text.x = element_text(angle=90,hjust=1,vjust = .5,face = 'plain',size=8))
ggsave("./output/figs/Fungal_Class_by_Size_Fraction_barplot.png",dpi=300,height = 8,width = 8)
x <- 
full_spp %>% 
  subset_taxa(Kingdom == "Fungi") %>% 
  transform_sample_counts(ra)

adonis2(formula = (x@otu_table %>% as('matrix')) ~ x@sam_data$size_fraction ) %>% 
  broom::tidy() %>% 
  saveRDS("./output/PermANOVA_Fungal_Community_vs_SizeClass.RDS")

plot_fungal_phylum <- 
function(x){
  x %>% 
    subset_taxa(Kingdom == "Fungi" & !is.na(Phylum)) %>% 
    transform_sample_counts(ra) %>% 
    plot_bar2(fill="Phylum") +
    scale_fill_viridis_d(option='turbo') +
    # scale_x_discrete(limits=depth_order) +
    facet_wrap(~size_fraction_um,scales='free')
}

p1 <- plot_fungal_phylum(ps_its) + ggtitle("ITS") + theme(axis.text.x = element_blank(), plot.title = element_text(color="blue",face='bold'))
p2 <- plot_fungal_phylum(ps_ssu) + ggtitle("SSU") + theme(axis.text.x = element_blank(), plot.title = element_text(color="blue",face='bold'))    
p3 <- plot_fungal_phylum(ps_lsu) + ggtitle("LSU")   + theme(axis.text.x = element_blank(), plot.title = element_text(color="blue",face='bold'))


full_spp@tax_table[,2] %>% unique %>% unname

p4 <- 
full_spp %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum)) %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Phylum") +
  scale_fill_viridis_d(option='turbo') +
  # scale_x_discrete(limits=depth_order) +
  facet_wrap(~size_fraction_um,scales='free') +
  ggtitle("Combined amplicons") +
  theme(plot.title = element_text(color="blue",face='bold'))

p1 / p2 / p3 / p4 + plot_layout(guides = 'collect')
ggsave("./output/figs/fungal_phylum_plots.png",height = 12,width = 8)


# how different do communities look between amplicons?


# FUNGI ONLY ####

fungi@tax_table[,1]




# glom taxa at species level?
# merge together?
# convert to relabund?

# METADATA PLOTS ####
ps_its@sam_data %>% names


# TAXONOMY PLOTS ####
ps_its %>% 
  merge_samples("size_fraction_um",fun = sum) %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Kingdom") +
  scale_fill_viridis_d(option = 'turbo')  + ggtitle("ITS1")
ggsave("./output/figs/ITS_kingdom_plot.png",height = 8, width = 6, dpi = 300)

ps_ssu %>% 
  merge_samples("size_fraction_um",fun = sum) %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Kingdom") +
  scale_fill_viridis_d(option = 'turbo') + ggtitle("SSU")
ggsave("./output/figs/SSU_kingdom_plot.png",height = 8, width = 6, dpi = 300)


ps_lsu %>% 
  merge_samples("size_fraction_um",fun = sum) %>% 
  transform_sample_counts(ra) %>% 
  plot_bar2(fill="Kingdom") +
  scale_fill_viridis_d(option = 'turbo')  + ggtitle("LSU")
ggsave("./output/figs/LSU_kingdom_plot.png",height = 8, width = 6, dpi = 300)

its <- 
  ps_its %>% 
  merge_samples("size_fraction_um",fun = sum) %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% mutate(amplicon = "ITS1")
ssu <- 
  ps_ssu %>% 
  merge_samples("size_fraction_um",fun = sum) %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% mutate(amplicon = "SSU")
lsu <- 
  ps_lsu %>% 
  merge_samples("size_fraction_um",fun = sum) %>% 
  transform_sample_counts(ra) %>% 
  psmelt() %>% mutate(amplicon = "LSU")

its %>% 
  full_join(ssu) %>% 
  full_join(lsu) %>% 
  ggplot(aes(x=Sample,y=Abundance,fill=Kingdom)) +
  geom_col() +
  facet_wrap(~amplicon) +
  scale_fill_viridis_d(option = 'turbo') +
  theme(legend.title = element_text(face='bold',size=14),
        legend.text = element_text(face='bold',size=12),
        axis.text.x = element_text(angle=90,hjust=1,vjust = .5)) +
  labs(x="Size fraction",y="Relative abundance")
ggsave("./output/figs/Kingdoms_by_amplicon_and_sizeclass.png",height = 8,width = 8,dpi=300)
