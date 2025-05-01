# 1. Preparation ----
## 1.1. Set WD ----
setwd("~/Projects/2023003_22057_PetrosiidaeMicrobiome/")

## 1.2. Package loading ----
library(phyloseq); packageVersion('phyloseq')
library(microbiome); packageVersion('microbiome')
library(DT)
library(vegan); packageVersion('vegan')
library(data.table)
library(pairwiseAdonis)
library(ape)
library(phytools)
library(dendextend)
library(phangorn)
library(ggtree)
library(rcompanion)
library(multcompView)
# library(picante)
# library(RVAideMemoire)
# library(TreeDist)
# library(ggplot2)
# library(ggdendro)
# library(plotly)
# library(plyr)
# library(reshape2)
# library(microViz)

## 1.3. Import data ---- 
ps = readRDS("ps_phylosymbiosis.rds")
ps

## 1.4 Merge samples by species ----
ps_merge = merge_samples(ps, "species")
ps_merge

## 1.5 Make compositional and rarefied datasets ----
ps_compo <- transform(ps, "compositional")   
ps_compo

sum(sample_sums(ps))
summary(sample_sums(ps))
ps_rarefied = rarefy_even_depth(ps)
ps_rarefied

## 1.6 visualize dataframe ----
metadata = as(sample_data(ps_compo), "data.frame")
datatable(metadata)

## 1.7 Set color vectors ----
color_species = c(
  "Acanthostrongylophora ingens" = "#ff8c00",
  "Neopetrosia carbonaria" = "#3cb371",
  "Neopetrosia chaliniformis" = "#0000ff",
  "Neopetrosia eurystomata" = "#000080",
  "Neopetrosia ovata" = "#f5deb3",
  "Neopetrosia proxima" = "#ff00ff",
  "Neopetrosia rosariensis" = "#dda0dd",
  "Petrosia (Petrosia) aff. elephantotus" = "#ba55d3",
  "Petrosia (Petrosia) elephantotus" = "#808000",
  "Petrosia (Petrosia) lignosa" = "#008080",
  "Petrosia (Petrosia) nigricans" = "#fa8072",
  "Petrosia (Petrosia) nova spec Curacao" = "#ff0000",
  "Petrosia (Petrosia) nova spec Lanyu" = "#8b0000",
  "Petrosia (Petrosia) weinbergi" = "#ffd700",
  "Petrosia (Strongylophora) corticata" = "#1e90ff",
  "Xestospongia mamillata" = "#ff1493",
  "Xestospongia muta" = "#00ff7f",
  "Xestospongia nova spec Lanyu" = "#87cefa",
  "Xestospongia testudinaria" = "#00ffff",
  "Xestospongia vansoesti" = "#adff2f",
  "Xestospongia viridenigra" = "#000000"
)

color_families = c("Archaea | Crenarchaeota | Nitrososphaeria | Nitrosopumilales | Nitrosopumilaceae"="#003300",                
                   "Bacteria | Acidobacteriota | Acidobacteriae | PAUC26f | NA"="#38B6FF",                                            
                   "Bacteria | Acidobacteriota | Thermoanaerobaculia | Thermoanaerobaculales | Thermoanaerobaculaceae"="#00FFE8",     
                   "Bacteria | Acidobacteriota | Vicinamibacteria | Subgroup 9 | NA"="#5CE1E6",                                       
                   "Bacteria | Acidobacteriota | Vicinamibacteria | Vicinamibacterales | NA"="#38B6FF",                               
                   "Bacteria | Actinobacteriota | Acidimicrobiia | Microtrichales | Microtrichaceae"="#000033",                       
                   "Bacteria | Bacteroidota | Rhodothermia | Rhodothermales | Rhodothermaceae"="#99CC99",
                   "Bacteria | Chloroflexi | Anaerolineae | Caldilineales | Caldilineaceae"="#F60000",                                
                   "Bacteria | Chloroflexi | Anaerolineae | SBR1031 | A4b"="#C20F0F",                                                 
                   "Bacteria | Chloroflexi | Dehalococcoidia | SAR202 clade | NA"="#FD8181",                                          
                   "Bacteria | Chloroflexi | TK10 | NA | NA"="#F60000",                                                               
                   "Bacteria | Chloroflexi | TK17 | NA | NA"="#B40202",                                                               
                   "Bacteria | Cyanobacteria | Cyanobacteriia | Synechococcales | Cyanobiaceae"="#00FF00",                            
                   "Bacteria | Dadabacteria | Dadabacteriia | Dadabacteriales | NA"="#ff9933",                                        
                   "Bacteria | Entotheonellaeota | Entotheonellia | Entotheonellales | Entotheonellaceae"="darkslateblue",
                   "Bacteria | Gemmatimonadota | BD2-11 terrestrial group | NA | NA"="#666600",                                       
                   "Bacteria | Myxococcota | bacteriap25 | NA | NA"="#FF66C4",                                                        
                   "Bacteria | Nitrospirota | Nitrospiria | Nitrospirales | Nitrospiraceae"="#D9D9D9",                                
                   "Bacteria | PAUC34f | NA | NA | NA"="#FFF847",                                                                     
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Defluviicoccales | NA"="#660066",                              
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodobacterales | Rhodobacteraceae"="#990099",   
                   "Bacteria | Proteobacteria | Alphaproteobacteria | Rhodospirillales | Magnetospiraceae" = "purple",
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Burkholderiales | Burkholderiaceae"="#C1FF72",                             
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Ectothiorhodospirales | Ectothiorhodospiraceae"="#7ED957",      
                   "Bacteria | Proteobacteria | Gammaproteobacteria | EPR3968-O8a-Bc78 | NA"="#00BF63",                              
                   "Bacteria | Proteobacteria | Gammaproteobacteria | NA | NA"="#5EFFA8",                                           
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Nitrosococcales | Nitrosococcaceae"="#B5E8A6",                 
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | Endozoicomonadaceae"="#57FF48",
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Pseudomonadales | KI89A clade"="#99FF33",                       
                   "Bacteria | Proteobacteria | Gammaproteobacteria | Steroidobacterales | Woeseiaceae"="#00CC00",                    
                   "Bacteria | Proteobacteria | Gammaproteobacteria | UBA10353 marine group | NA"="#66FF00",                          
                   "Bacteria | Spirochaetota | Spirochaetia | Spirochaetales | Spirochaetaceae"="#724C1F",
                   "Bacteria | Verrucomicrobiota | Verrucomicrobiae | Pedosphaerales | Pedosphaeraceae" = "yellow",
                   "Other"="grey")

# 2. Analyses -------------------------------------------------------------------------
## 2.1 Rarefaction curves ----
Rarecurve = rarecurve(as.data.frame(t(otu_table(ps))), step = 100, col = "blue", label= FALSE)

#save manually the plot to ./1. Data_prep results

## 2.2 General microbiome dataset information ----

### Summarize dataset ----

summarize_phyloseq(ps)
ps_ASV = prune_taxa(taxa_sums(ps) > 0, ps) #Remove ASV's that are not present in the dataset
ps_ASV #Number of unique ASVs in dataset
  
### Microbial phyla ---- 
# To see general information about the microbial phyla present in the dataset and the amount of reads. Aggregate ASVs at the phylum level
ps_phylum <- tax_glom(ps_ASV, "Phylum")
ps_phylum
taxa_names(ps_phylum) <- tax_table(ps_phylum)[, 2]
taxa_names(ps_phylum)

phylum_reads = colSums(t(otu_table(ps_phylum)))
sort(phylum_reads, decreasing = TRUE) # Check the most abundant phyla and number of reads per phylum


### distribution of taxa ----
ps_df_taxa <- data.table(tax_table(ps_ASV),
                          ASVabundance = taxa_sums(ps_ASV),
                          ASV = taxa_names(ps_ASV))

ps_tax_plot <- ggplot(ps_df_taxa, aes(ASVabundance)) +
  geom_histogram() + ggtitle("Histogram of ASVs (unique sequence) counts") +
  theme_bw() + scale_x_log10() + ylab("Frequency of ASVs") + xlab("Abundance (raw counts)")

print(ps_tax_plot)

### reads per sample ----
ps_df= data.table(as(sample_data(ps_ASV), "data.frame"), Reads_per_sample = sample_sums(ps_ASV), keep.rownames = TRUE)
ps_df_plot = ggplot(ps_df, aes(Reads_per_sample)) + geom_histogram() + ggtitle("Distribution of reads per sample") + ylab("Sample counts") 

print(ps_df_plot) #normal plot

### Distribution of ASVs ----
ps.dt.taxa = data.table(tax_table(ps_ASV),OTUabundance = taxa_sums(ps_ASV),OTU = taxa_names(ps_ASV))
ps.dt.tax.plot <- ggplot(ps.dt.taxa, aes(OTUabundance)) + geom_histogram() + ggtitle("Histogram of OTU (unique sequence) counts") + theme_bw()
print(ps.dt.tax.plot)

### Variance ----
Variance.plot <- qplot(log10(apply(otu_table(ps), 1, var)), xlab = "log10(variance)", main = "Variance in OTUs")
print(Variance.plot)


## 2.3 Beta diversity ----
### 2.3.1 Make NMDS ----

# nmds_bray <- ordinate(ps_compo, "NMDS", "bray")
# nmds_bray$stress

# write.csv(nmds_bray$points , file.path("./3_Beta_div_results/","nmds_bray-points.csv"))
# saveRDS(nmds, file.path("./3_Beta_div_results/","nmds_bray.RDS"))
nmds_bray = readRDS(file.path("./3_Beta_div_results/","nmds_bray.RDS"))
nmds_bray$stress

### 2.3.2 Plot NMDS ----
data_nmds <- plot_ordination(ps_compo, nmds_bray, type = "Samples", justDF = TRUE)

plot_nmds = ggplot(data_nmds, aes(x = NMDS1, y = NMDS2))
plot_nmds = plot_nmds + geom_point(aes(fill = species.full), pch = 21, size = 3.5, alpha = 0.8)
# plot_nmds = plot_nmds + scale_shape_manual(values = shape_genus)
plot_nmds = plot_nmds + scale_fill_manual(values = color_species)
plot_nmds = plot_nmds + theme_bw(base_size = 15)+ theme(plot.title = element_text(size=22, face="bold"), axis.text=element_text(size=12),axis.title=element_text(size=14))
plot_nmds = plot_nmds + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7), ncol = 1))
plot_nmds = plot_nmds + annotate("text", label="2D Stress = 0.201", x = -0.25, y = 0.8)
# plot_nmds = plot_nmds + ggtitle("NMDS based on Bray-Curtis distances of sponge species")
plot_nmds = plot_nmds + labs(fill = "Species") + theme(legend.title = element_text(size = 17))
plot_nmds

ggsave(filename = "Plot_NMDS_bray.pdf",
       plot = plot_nmds, 
       device = "pdf" , 
       width = 35 , height = 20, units = "cm", 
       path = "./3_Beta_div_results")

### 2.3.3 PERMANOVA ----

permanova.result = adonis2(distance(ps_compo, method="bray") ~ species.full, data = metadata, by = NULL)
datatable(permanova.result)

write.csv(as.data.frame(permanova.result), file.path("./3_Beta_div_results/","permanova-result-species.csv"))

data_pairwiseadonis_species = pairwise.adonis(distance(ps_compo, method="bray"), metadata$species.full)
datatable(data_pairwiseadonis_species)

write.csv(as.data.frame(data_pairwiseadonis_species), file.path("./3_Beta_div_results/","data-pairwiseadonis-species.csv"))

## 2.4 Phylosymbiosis ----

### 2.4.1 Mantel ----
#For mantel test, the distance matrixes need to be sorted as they are expected to have the same layout. Exported Geneious distance matrixes are sorted by phylogenetic order and distance matrix from ps.object are sorted numerically. Make sure to edit phylo dist.matrix to sort the first column, then transpose and then sort again. This results in a usable genetic distance matrix sorted by alphabetical order, as required for a mantel test here.

#generate distance matrix microbial communities
dist.matrix.bray = as.matrix(distance(ps_compo, method = "bray"))
dist.matrix.bray
dim(dist.matrix.bray)

dist.matrix.wunifrac = as.matrix(distance(ps_compo, method = "wunifrac"))
dist.matrix.wunifrac
dim(dist.matrix.wunifrac)

dist.matrix.uunifrac = as.matrix(distance(ps_compo, method = "uunifrac"))
dist.matrix.uunifrac
dim(dist.matrix.uunifrac)

#import phylogenetic distance matrix
phylo.dist.matrix = read.csv(file = "./data/Barcoding data/28S RAxML Tree_distmatrix_sorted.csv" , header = TRUE , sep = "," , row.names = 1)
phylo.dist.matrix
dim(phylo.dist.matrix)

data_mantel_bray = mantel(dist.matrix.bray, phylo.dist.matrix, method = "spearman", permutations = 9999, na.rm = FALSE)
data_mantel_bray

data_mantel_wunifrac = mantel(dist.matrix.wunifrac, phylo.dist.matrix, method = "spearman", permutations = 9999, na.rm = FALSE)
data_mantel_wunifrac

data_mantel_uunifrac = mantel(dist.matrix.uunifrac, phylo.dist.matrix, method = "spearman", permutations = 9999, na.rm = FALSE)
data_mantel_uunifrac

### 2.4.2 RF & Entanglement ----
#make Dendrogram of community
micro.dend.bray <- hclust(phyloseq::distance(ps_compo, method="bray"))
micro.dend.bray

micro.dend.wunifrac <- hclust(phyloseq::distance(ps_compo, method="wunifrac"))
micro.dend.wunifrac

micro.dend.uunifrac <- hclust(phyloseq::distance(ps_compo, method="uunifrac"))
micro.dend.uunifrac

#import phylogenetic tree
phylo.tree = read.tree("./data/Barcoding data/28S RAxML Tree.newick")
phylo.tree$tip.label #check if tip labels match the sample numbers

plot(phylo.tree) #check tree

#Calculate nRF values.
nRF.bray = RF.dist(phylo.tree,as.phylo(micro.dend.bray),normalize = TRUE)
nRF.bray

nRF.wunifrac = RF.dist(phylo.tree,as.phylo(micro.dend.wunifrac),normalize = TRUE)
nRF.wunifrac

nRF.uunifrac = RF.dist(phylo.tree,as.phylo(micro.dend.uunifrac),normalize = TRUE)
nRF.uunifrac

#Calculate RF with p-values

RF.bray = cospeciation(phylo.tree, as.phylo(micro.dend.bray))
RF.bray

RF.wunifrac = cospeciation(phylo.tree, as.phylo(micro.dend.wunifrac))
RF.wunifrac

RF.uunifrac = cospeciation(phylo.tree, as.phylo(micro.dend.uunifrac))
RF.uunifrac


### 2.4.3 Tanglegram plot ----
#(Optional) Update labels for dendrogram plot
# micro.dend$labels
# micro.dend$labels <- paste(metadata$extract.ID, metadata$species, metadata$fieldnr.suffix)[match(micro.dend$labels, rownames(metadata))]
# 
# tree$tip.label
# tree$tip.label <- paste(metadata$extract.ID, metadata$species, metadata$fieldnr.suffix)[match(tree$tip.label, rownames(metadata))]

#prepare tree and dendrogram for tanglegram
phylo.tree.root = midpoint_root(phylo.tree)
phylo.tree.root.ultra = force.ultrametric(phylo.tree.root, method=c("nnls"))
plot(phylo.tree.root.ultra) #check if tree is ultrametric

dends.bray = dendlist(phylo.tree.root.ultra, micro.dend.bray)
dends.bray.untangled = untangle(dends.bray, method = "step2side")
entanglement.bray = entanglement(dends.bray.untangled)
entanglement.bray

dends.wunifrac = dendlist(phylo.tree.root.ultra, micro.dend.wunifrac)
dends.wunifrac.untangled = untangle(dends.wunifrac, method = "step2side")
entanglement.wunifrac = entanglement(dends.wunifrac.untangled)
entanglement.wunifrac

dends.uunifrac = dendlist(phylo.tree.root.ultra, micro.dend.uunifrac)
dends.uunifrac.untangled = untangle(dends.uunifrac, method = "step2side")
entanglement.uunifrac = entanglement(dends.uunifrac.untangled)
entanglement.uunifrac

tanglegram.plot = tanglegram(dends.bray.untangled, 
               highlight_distinct_edges = FALSE, 
               common_subtrees_color_lines = FALSE, 
               highlight_branches_lwd = FALSE, 
               lwd = 1,
               main_left = "Sponge phylogeny",
               main_right = "Microbial dendrogram",
               main = "28S_raxml_bray",
               cex_main = 1,
               columns_width = c(2,2,2),
               # left_dendo_mar = c(1,1,1,3),
               # right_dendo_mar = c(1,3,1,1),
               # margin_inner = 3,
               # margin_outer = 3,
               lab.cex = 0.5,
               dLeaf_right = FALSE,
)
tanglegram.plot

#save manually the plot to ./5_phylosymbiosis_results. Manually edited in Adobe illustrator for publishable graph. Added depth and geography metadata and coloured based on species.

### 2.4.4 Phylogenetic tree and phylosignal ----

# Calculate Shannon diversity for each sampple
shannon<-vegan::diversity(t(otu_table(ps_rarefied)), index = "shannon") #for merged dataset, no tranfsormation of OTU table
shannon

#Calculate Phylogenetic signal on Alpha diversity for every sample
pagels<-phylosig(phylo.tree,shannon,method="lambda",test=T)
pagels


# Import consensus tree

phylo.tree.consensus = read.tree("./data/Barcoding data/28S_consensus RAxML Tree.newick")
phylo.tree.consensus
plot(phylo.tree.consensus) #Check tree
phylo.tree.consensus$tip.label
phylo.tree.consensus$tip.label <- gsub(" _consensus", "", phylo.tree.consensus$tip.label) #for consensus tree
phylo.tree.consensus$tip.label <- gsub("'", "", phylo.tree.consensus$tip.label) #for refseq tree named
phylo.tree.consensus$tip.label
plot(phylo.tree.consensus) #Check tree

phylo.tree.consensus.root = midpoint_root(phylo.tree.consensus)
phylo.tree.consensus.root.ultra = force.ultrametric(phylo.tree.consensus.root, method=c("nnls"))
plot(phylo.tree.consensus.root.ultra) #check if tree is ultrametric


#Calculate average shannon diversity per species
metadata$shannon <- shannon[match(rownames(metadata), names(shannon))] #add shannon values to metadata column

metadata.shannon.consensus = aggregate(shannon ~ species, data = metadata, FUN = mean, na.rm = TRUE)
datatable(metadata.shannon.consensus)
metadata.shannon.consensus
dim(metadata.shannon.consensus)


#Calculate Phylogenetic signal on Alpha diversity for average shannon diversity per species
shannon.consensus <- setNames(metadata.shannon.consensus$shannon, metadata.shannon.consensus$species)
shannon.consensus

#verify tip labels of phylo.tree.consensus (for whatever reasons tip labels and speices names don't match, even if they do)
phylo.tree.consensus$tip.label
phylo.tree.consensus$tip.label[9] = "Neopetrosia rosariensis"
names(shannon.consensus)[names(shannon.consensus) == "Neopetrosia rosariensis"] = "Neopetrosia rosariensis"

setdiff(phylo.tree.consensus$tip.label, names(shannon.consensus))
setdiff(names(shannon.consensus), phylo.tree.consensus$tip.label)

pagels.consensus<-phylosig(phylo.tree.consensus,shannon.consensus,method="lambda",test=T)
pagels.consensus

#Plot consensus tree for Figure 1

# ggtree(phylo.tree.consensus.root.ultra) +
#   geom_tippoint(aes(size = metadata.shannon.consensus$shannon), color = "blue", alpha = 0.7) +
#   scale_size(range = c(2, 10)) +  # Adjust circle sizes
#   theme_minimal() +
#   labs(size = "Shannon Diversity")  # Legend label

consensus.tree.plot = ggtree(phylo.tree.consensus.root.ultra)
consensus.tree.plot = consensus.tree.plot + geom_tiplab()
consensus.tree.plot = consensus.tree.plot + geom_text2(aes(subset = !isTip, label = label), hjust = -0.3)
consensus.tree.plot = consensus.tree.plot + xlim_tree(1.5)
consensus.tree.plot

datatable(consensus.tree.plot$data)

consensus.tree.plot$data

consensus.tree.plot$data <- merge(consensus.tree.plot$data, 
                                  metadata.shannon.consensus, 
                                  by.x = "label",  # Match species names in tree
                                  by.y = "species", 
                                  all.x = TRUE)  # Keep all tree tips, even if Shannon is missing

# Plot tree with Shannon diversity visualized as tip point size
consensus.tree.plot.shannon = consensus.tree.plot + 
  geom_tippoint(aes(size = shannon), color = "blue", alpha = 0.7) +  
  scale_size(range = c(3, 10), name = "Shannon Diversity") +  
  theme_void()
consensus.tree.plot.shannon

ggsave(filename = "Plot_consensus.tree.shannon.pdf",
       plot = consensus.tree.plot.shannon, 
       device = "pdf" , 
       width = 20 , height = 20, units = "cm", 
       path = "./2_Alpha_div_results")

## For editing in Adobe Illustrator for figure 1. Combined with Community Barplot

### 2.4.5 Community barplot ----

ps_merge
sample_names(ps_merge)

ps_merge_compo <- transform(ps_merge, "compositional")   
ps_merge_compo

plot_bar_2 <-  function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, facet_grid = NULL, border_color = NA) 
{
  mdf = psmelt(physeq)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack",  color = border_color)
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


# Define your custom order of samples (character string of sample names)
phylo.consensus.order <- c(
  "Neopetrosia ovata",
  "Neopetrosia rosariensis",
  "Petrosia Petrosia weinbergi",
  "Petrosia Petrosia nova spec Lanyu",
  "Petrosia Petrosia nova spec Curacao",
  "Petrosia Strongylophora corticata",
  "Petrosia Petrosia elephantotus",
  "Petrosia Petrosia aff elephantotus",
  "Xestospongia viridenigra",
  "Neopetrosia carbonaria",
  "Xestospongia nova spec Lanyu",
  "Xestospongia muta",
  "Xestospongia testudinaria",
  "Neopetrosia eurystomata",
  "Neopetrosia proxima",
  "Neopetrosia chaliniformis",
  "Petrosia Petrosia lignosa",
  "Petrosia Petrosia nigricans",
  "Xestospongia mamillata",
  "Acanthostrongylophora ingens",
  "Xestospongia vansoesti"
)

physeq_aggreg <- aggregate_rare(ps_merge_compo, level = "Family", detection = 3/100, prevalence = 0/100)
physeq_aggreg

ntaxa(physeq_aggreg)
tax.table.family = as.data.frame(tax_table(physeq_aggreg))
tax.table.family$Family


sample_data(physeq_aggreg)$consensus.order <- factor(sample_names(physeq_aggreg))
levels(sample_data(physeq_aggreg)$consensus.order)
sample_data(physeq_aggreg)$consensus.order <- factor(sample_data(physeq_aggreg)$consensus.order, levels = rev(phylo.consensus.order))

barplot_family = plot_bar_2(physeq_aggreg, "consensus.order", fill = "Family")
barplot_family = barplot_family + scale_fill_manual(values = color_families)
barplot_family = barplot_family +  coord_flip()
barplot_family = barplot_family + theme(legend.text = element_text(size=9),
                                        # axis.text.x = element_blank(),
                                        legend.title = element_blank(),
                                        axis.ticks.x=element_blank(), 
                                        axis.title = element_blank())
barplot_family = barplot_family + guides (fill = guide_legend(ncol = 1)) 
barplot_family
barplot_family = barplot_family + theme(panel.spacing = unit(0, "cm", data = NULL),panel.border = element_rect(color = "black", fill = NA, size = 0.9))
barplot_family = barplot_family +scale_y_continuous(expand = c(0,0))
barplot_family

ggsave(filename = "barplot_28S.consensus_orderphylogeny.pdf", 
       plot = barplot_family, 
       device = "pdf" , 
       width = 40 , height = 21 , units = "cm", 
       path = "./4_Compositional_results")

#merge consensus phylogeny with average shannon value and barplot for figure 1

## 2.5 Alpha Diversity ----
data_alpha = estimate_richness(ps_rarefied , measures = c("Observed","chao1", "Shannon"))
data_alpha

Pielou = data_alpha$Shannon / log(data_alpha$Observed)
Pielou

data_alpha_all = cbind(metadata["species.full"], data_alpha , Pielou)
data_alpha_all

write.csv(data_alpha_all, file.path("./2_Alpha_div_results" , "data_alpha_all.csv"))

data_long = reshape2::melt(data_alpha_all, id.var= "species.full")
data_long

#Plot species
plot_alpha_all = ggplot(data_long, aes(x = species.full, y = value))
plot_alpha_all = plot_alpha_all + geom_point(aes(fill = species.full), size = 2, alpha = 0.8, pch = 21)
plot_alpha_all = plot_alpha_all + geom_boxplot(aes(fill = species.full), alpha = 0.8, size = 0.8)
plot_alpha_all = plot_alpha_all + scale_fill_manual(values = color_species)
plot_alpha_all = plot_alpha_all + theme_bw(base_size = 15)  + theme(legend.position="left")
plot_alpha_all = plot_alpha_all + facet_grid(variable ~ . , scales = "free")
plot_alpha_all = plot_alpha_all + scale_shape_manual(values = c(21,22,23,24,25))
plot_alpha_all = plot_alpha_all + theme(axis.title.x = element_blank(),axis.title.y = element_blank())
plot_alpha_all = plot_alpha_all + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1))
plot_alpha_all = plot_alpha_all + guides(fill = guide_legend(override.aes = list(shape = 21, size = 7), ncol = 1))
plot_alpha_all = plot_alpha_all + labs(fill = "Species") + theme(legend.title = element_text(size = 17))
plot_alpha_all

ggsave(filename = "Alpha diversity.pdf",
       path = "~/Projects/2023003_22057_PetrosiidaeMicrobiome/2_Alpha_div_results",
       plot = plot_alpha_all, 
       width = 12,
       height = 10)

## Test for normality of the distribution of the indexes
# P values <0.05 imply no normal distribution --> thus Kruskall-wallis.
# If P value >0.05 use Parametric ANOVA
shapiro_data_shannon = shapiro.test(data_alpha_all$Shannon)
shapiro_data_shannon

shapiro_data_chao1 = shapiro.test(data_alpha_all$Chao1)
shapiro_data_chao1

shapiro_data_pielou = shapiro.test(data_alpha_all$Pielou)
shapiro_data_pielou

shapiro_data_alpha <- matrix(nrow = 3 ,  ncol=2, byrow=TRUE)
colnames(shapiro_data_alpha) = c("W","p-value")
rownames(shapiro_data_alpha) = c("Shannon","Chao1","Pielou")

shapiro_data_alpha[,1] <- c(shapiro_data_shannon$statistic,
                            shapiro_data_chao1$statistic,
                            shapiro_data_pielou$statistic)

shapiro_data_alpha[,2]<- c(shapiro_data_shannon$p.value,
                           shapiro_data_chao1$p.value,
                           shapiro_data_pielou$p.value)

datatable(shapiro_data_alpha)

write.csv(shapiro_data_alpha, file.path("./2_Alpha_div_results" , "Shapiro_data_alpha.csv"))

#Test for differences between species 
data_kruskal_Shannon_species = kruskal.test(Shannon ~ species.full, data_alpha_all)
data_kruskal_Shannon_species

data_kruskal_Chao1_species = kruskal.test(Chao1 ~ species.full, data_alpha_all)
data_kruskal_Chao1_species

data_kruskal_Pielou_species = kruskal.test(Pielou ~ species.full, data_alpha_all)
data_kruskal_Pielou_species 

data_kruskal_alpha_species  <- matrix(nrow = 3 ,  ncol=3, byrow=TRUE)
colnames(data_kruskal_alpha_species) = c("Chi-square","Df","p-value")
rownames(data_kruskal_alpha_species) = c("Shannon","Chao1","Pielou")

data_kruskal_alpha_species

data_kruskal_alpha_species[,1] <- c(data_kruskal_Shannon_species$statistic,
                                  data_kruskal_Chao1_species$statistic,
                                  data_kruskal_Pielou_species$statistic)

data_kruskal_alpha_species[,2]<- c(data_kruskal_Shannon_species$parameter,
                                 data_kruskal_Chao1_species$parameter,
                                 data_kruskal_Pielou_species$parameter)

data_kruskal_alpha_species[,3]<- c(data_kruskal_Shannon_species$p.value,
                                 data_kruskal_Chao1_species$p.value,
                                 data_kruskal_Pielou_species$p.value)

data_kruskal_alpha_species

write.csv(data_kruskal_alpha_species, file.path("./2_Alpha_div_results" , "Data_kruskal_alpha_species.csv"))

#After kruskal wallis - do pairwise test with Wilcoxon 
data_wilcox_shannon_species = pairwise.wilcox.test(data_alpha_all$Shannon, data_alpha_all$species.full, p.adjust.method = "bonferroni")
data_wilcox_shannon_species
data_wilcox_shannon_species$p.value



data_wilcox_shannon_species_full = fullPTable(data_wilcox_shannon_species$p.value)
data_wilcox_shannon_species_full



indices_wilcox_shannon_species = multcompLetters(data_wilcox_shannon_species_full,
                                           compare="<",
                                           threshold=0.05,
                                           Letters=letters,
                                           reversed = FALSE)
indices_wilcox_shannon_species$Letters



data_wilcox_shannon_species_full2 = cbind(data_wilcox_shannon_species_full, indices_wilcox_shannon_species$Letters )
data_wilcox_shannon_species_full2

datatable(data_wilcox_shannon_species_full2)

write.table(data_wilcox_shannon_species_full2, 
            file.path("./2_Alpha_div_results", "Data_wilcox_shannon_species.csv"), 
            sep = ";", 
            dec = ",", 
            row.names = TRUE, 
            col.names = NA)


## 2.6 Compare phylogenies ----

#Export trees as newick and import here. Samples should have the same name (enumber)

### 2.6.1 Compare 28S to Concatenated ----
#Import trees to compare
treeconcat = read_tree("./data/Barcoding data/concatenated_enum RAxML Tree.newick") #all concatenated sequences
treeconcat.28S = read_tree("./data/Barcoding data/28S_concatenated_enum RAxML Tree.newick") #28S tree using only samples used in concatenated tree

treeconcat$tip.label #check if tip labels are matching
treeconcat.28S$tip.label #check if tip labels are matching

## Dendextend (Only works when trees have same amount of tips)

#prepare tree and dendrogram for tanglegram
treeconcat.28S.root = midpoint_root(treeconcat.28S)
treeconcat.28S.root.ultra = force.ultrametric(treeconcat.28S.root, method=c("nnls"))
plot(treeconcat.28S.root.ultra) #check if tree is ultrametric

treeconcat.root = midpoint_root(treeconcat)
treeconcat.root.ultra = force.ultrametric(treeconcat.root, method=c("nnls"))
plot(treeconcat.root.ultra) #check if tree is ultrametric

dends.phylo = dendlist(treeconcat.root.ultra, treeconcat.28S.root.ultra)
dends.phylo.untangled = untangle(dends.phylo, method = "step2side")
entanglement.phylo = entanglement(dends.phylo.untangled)
entanglement.phylo

tanglegram.phylo.plot = tanglegram(dends.phylo.untangled, 
                             highlight_distinct_edges = FALSE, 
                             common_subtrees_color_lines = FALSE, 
                             highlight_branches_lwd = FALSE, 
                             lwd = 1,
                             main_left = "Treeconcat",
                             main_right = "Treeconcat.28S",
                             main = "tanglegram concatenated vs 28S",
                             cex_main = 1,
                             columns_width = c(2,2,2),
                             # left_dendo_mar = c(1,1,1,3),
                             # right_dendo_mar = c(1,3,1,1),
                             # margin_inner = 3,
                             # margin_outer = 3,
                             lab.cex = 0.5,
                             dLeaf_right = FALSE,
)
tanglegram.phylo.plot

#Save manually to  "tanglegram_28S.concatenated28S.pdf" For supplementary figure S2
