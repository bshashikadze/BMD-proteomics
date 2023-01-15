myocardium
================
BS
2023-01-15

DMD -\> “Δ52” 
BMD -\> “Δ51-52”

“proteinGroups_myocardium.txt” can be downloaded from PRIDE repository
(ID to be provided); all other necessary files can be found in this repository

## load libraries

``` r
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(ggforce)
library(cowplot)
library(grid)
library(circlize)
library(WebGestaltR)
library(xlsx)
library(msigdbr)
```

## load protein groups file (MaxQuant output)

``` r
proteingroups <- read.delim("proteinGroups_myocardium.txt", sep = "\t", header = T)
```

## data filtering

### filter protein groups and write conditions file (Perseus -\> filter rows based on categorical column)

only identified by site, potential contaminants and reverse hits will be
removed modify conditions file according to the experimental groups

``` r
filter_pg_fn <- function(data) {
  
  # filter
  data <- data %>% 
        filter(Potential.contaminant   != "+") %>%
        filter(Only.identified.by.site != "+") %>%
        filter(Reverse != "+") %>% 
        select(Protein.IDs, starts_with("LFQ.Intensity."))
 
  # create conditions file    
 Bioreplicate <- colnames(data)[! colnames(data) %in% c('Protein.IDs')]
 Group    <- str_remove(Bioreplicate, "LFQ.intensity.")
 Conditions   <- data.frame(Bioreplicate, Group) %>%
 write.table("conditions.txt", row.names = F, sep = "\t")
 cat("conditions file was generated, rename to conditions_modified and modify the second column
    according to experimental groups")
 return(data)
}
# apply function
PG_filtered <- filter_pg_fn(proteingroups)
```


``` r
conditions <- read.delim("conditions_modified.txt", sep = "\t", header = T)
```

### filter for valid values (Perseus -\> filter rows based on valid values)

``` r
filter_valid_fn <- function(data, conditions_file, subset_groups = c(), 
                          min_val_atleastonegr =3) {
  
  # subset 
  data_subset <- data %>% 
    pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -Protein.IDs) %>% 
    left_join(conditions_file) %>% 
    filter(Group %in% subset_groups)
  
  
  # count conditions and subset (if necessary)
  n_valid <- data_subset %>% 
    group_by(Protein.IDs, Group) %>% 
    summarise(n_cond = sum(Intensity > 0, na.rm = TRUE)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = "Group", values_from = "n_cond", Protein.IDs)
  
  # filter for valid values
  data_subset_filtered <- data_subset %>% 
      pivot_wider(names_from = "Bioreplicate", values_from = "Intensity", Protein.IDs) %>% 
      left_join(n_valid) %>% 
      filter_at(vars(-starts_with("LFQ."), -starts_with("Protein.IDs")), 
                any_vars(. >= min_val_atleastonegr)) %>% 
      select(Protein.IDs, starts_with("LFQ.")) %>% 
      column_to_rownames("Protein.IDs") %>% 
      mutate_all(~na_if(., 0)) %>% 
      mutate_all(., log2) 
  
  return(data_subset_filtered)}
# apply function
PG_valid <- filter_valid_fn(PG_filtered, conditions_file = conditions, 
                            subset_groups = c("Δ51-52", "Δ52", "WT"), 
                            min_val_atleastonegr = 3)
```

## impute missing values from normal distribution (Perseus -\> impute missing values)

``` r
impute_fn <- function(data, downshift = 1.8, width = 0.3) {
  
   # obtain statistics from the data to build a distribution from which values will be imputed
    valid_data_descriptives <- data %>%
    rownames_to_column("id") %>%
    pivot_longer(names_to  = "Bioreplicate", values_to = "Intensity", -id) %>%
    group_by(Bioreplicate) %>%
    summarise(mean_valid   = mean(Intensity, na.rm = T),
              median_valid = median(Intensity, na.rm = T),
              sd_valid     = sd(Intensity,   na.rm = T),
              n_valid      = sum(!is.na(Intensity)),
              n_missing    = nrow(data) - n_valid) %>%
    ungroup()
    
    
  # impute missing values
  # imputation
  column_indices <- 1:nrow(valid_data_descriptives)
  random_data <- list()
  
  # makes the list which contains as many elements as the samples are and as many random values as the     
  # missing values are in each sample
  for (i in column_indices) {
    set.seed(seed = 1234)
    random_data[[i]] <- rnorm(n = valid_data_descriptives$n_missing[i],
                              mean = valid_data_descriptives$median_valid[i] - (downshift * valid_data_descriptives$sd_valid[i]),
                              sd = valid_data_descriptives$sd_valid[i]   * width)
    
  # impute the missing values
    data[is.na(data)[, valid_data_descriptives$Bioreplicate[i]],
              valid_data_descriptives$Bioreplicate[i]] <- random_data[[i]]}
    return(data)
}
# apply function
PG_imputed <- impute_fn(PG_valid)
```

## save data for statistical analysis in perseus

perform 1 way anova followed by permutation based FDR correction in
Perseus (FDR\<0.05;S0=0) perform t-test followed by permutation based
FDR correction in Perseus (FDR\<0.05;S0=0.1) Perseus version 1.6.15.0

note: missing value imputation is performed only once on original data

``` r
# for 1 way anova
write.table(PG_imputed %>% rownames_to_column("Protein_IDs"), "one_way_anova.txt", sep = "\t", row.names = F, quote = F)
# DMD versus WT
DMD_WT <- filter_valid_fn(PG_filtered, conditions_file = conditions, subset_groups = c("Δ52", "WT"), min_val_atleastonegr = 3) %>% 
  select(-starts_with("LFQ")) %>% 
  rownames_to_column("Protein_IDs") %>% 
  left_join(PG_imputed %>% 
              select(starts_with(conditions$Bioreplicate[conditions$Group != "Δ51-52"])) %>% 
              rownames_to_column("Protein_IDs"))

# save
write.table(DMD_WT, "DMD_WT.txt", sep = "\t", row.names = F, quote = F)
# BMD versus WT
BMD_WT <- filter_valid_fn(PG_filtered, conditions_file = conditions, subset_groups = c("Δ51-52", "WT"), min_val_atleastonegr = 3) %>% 
  select(-starts_with("LFQ")) %>% 
  rownames_to_column("Protein_IDs") %>% 
  left_join(PG_imputed %>% 
              select(starts_with(conditions$Bioreplicate[conditions$Group != "Δ52"])) %>% 
              rownames_to_column("Protein_IDs"))

write.table(BMD_WT, "BMD_WT.txt", sep = "\t", row.names = F, quote = F)
```

## read outputs of Perseus

### 1 way anova results

``` r
# 1 way anova results
anova_sig <- read.delim("one_way_anova_perseus.txt", sep = "\t", header = T) %>% 
  filter(C..ANOVA.Significant == "+") %>% 
  column_to_rownames("T..Protein_IDs") %>% 
  select(starts_with("LFQ")) %>% 
  rename_all(~str_replace(., "LFQ.intensity.", ""))
# standardization
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}
 
anova_sig <- t(apply(anova_sig, 1, cal_z_score))
```

### results of the t-test with x and y coordinates of s0 curve

``` r
# DMD vs WT
DMD_WT_perseus <- read.delim("DMD_WT_perseus.txt", sep = "\t", header = T)
xy_DMD_WT      <- read.delim("xy_DMD_WT.txt", sep = "\t", header = T)
# BMD vs WT
BMD_WT_perseus <- read.delim("BMD_WT_perseus.txt", sep = "\t", header = T)
xy_BMD_WT      <- read.delim("xy_BMD_WT.txt", sep = "\t", header = T)
```

## Hierarchical clustering of anova significant results

### plot heatmap

``` r
#set colors
col_fun = colorRamp2(c(-1, 0, 1), c("#56B4E9", "white",  "firebrick3"))
colours <- list('Group' = c('Δ52' = "#708090", 'WT' = "navy", 'Δ51-52' = "#D55E00"))

# column annotation
colAnn <- HeatmapAnnotation(df = conditions %>% select(Group), show_legend = F, 
                            which = 'col',
                            annotation_name_gp= gpar(fontsize = 9),
                            annotation_name_side = "left",
                            col = colours,
                            height = unit(0.4, "cm"),
                            simple_anno_size_adjust = TRUE)

# row annotation
rowAnn <- rowAnnotation(foo = anno_block(gp = gpar(fill = 'azure2'),
                              labels = c("Cl 1", "Cl 2"), 
                              width = unit(0.3, "cm"),
                              labels_gp = gpar(col = "black", 
                              fill = "white", fontsize = 9)))


# k-means clustering
pa = cluster::pam(as.matrix(anova_sig) , k = 2)


# ploting
set.seed(1234)
Hmap         <- Heatmap(as.matrix(anova_sig) ,
                top_annotation=colAnn, 
                show_heatmap_legend = F,
                show_row_names = F, 
                show_column_names = T,
                left_annotation = rowAnn,
                clustering_distance_rows    = "euclidean",
                clustering_distance_columns = "euclidean",
                clustering_method_columns   = "average",
                clustering_method_rows      = "average",
                col = col_fun, 
                row_dend_width     = unit(0.4, "cm"),
                column_dend_height = unit(0.4, "cm"),
                row_title = NULL,
                row_split = pa$clustering,
                width              = unit(1.5, "in"),
                height             = unit(3.49, "in"), 
                border = TRUE, column_names_gp = grid::gpar(fontsize = 9),
                row_names_gp = grid::gpar(fontsize = 9))


#calculate actual plot size
ht <- draw(Hmap)


w1 = ComplexHeatmap:::width(ht)
w1 = convertX(w1, "inch", valueOnly = TRUE)
h1 = ComplexHeatmap:::height(ht)
h1 = convertY(h1, "inch", valueOnly = TRUE)
c(w1, h1)


# for cowplot
heatmap_plot = grid.grabExpr(draw(Hmap))
```

### get the heatmap legend (more convenient to arrange e.g. in inkscape)

``` r
# ploting legends separatelly
# intensity legend
hm_legend = grid.grabExpr(color_mapping_legend(Hmap@matrix_color_mapping, 
                                               plot = T, title = "Intensity", 
                                               title_position = "topcenter", 
                                               legend_direction = c("horizontal"),  
                                               title_gp = gpar(fontsize = 8.5), 
                                               param = list(at = c(-1,  1), 
                                                            labels = c("low", "high"), 
                                                            legend_width = unit(2, "cm"), 
                                                            labels_gp = gpar(fontsize=8.5)),  
                                               labels_gp = gpar(fontsize = 7)))
# annotation legend
Group = grid.grabExpr(color_mapping_legend(Hmap@top_annotation@anno_list[["Group"]]@color_mapping, 
                                           nrow = 3, plot = T,  
                                           title_gp = gpar(fontsize = 8.5),  
                                           labels_gp = gpar(fontsize = 8.5)))

# combine legends and save
plot_grid(hm_legend, Group, ncol = 3, rel_widths = c(1,1,1), rel_heights = c(1,1,1))

ggsave("legend_hm.svg", width =3, height = 1)
```

## assign cluster number to each entry, plot line plot and perform ORA

### assign cluster number to each protein (this step should be perfomed carefully as the k-means order is not always corresponding to mathematical order (1,2,3..), best way is to delete row_title = NULL in complexheatmap and look the original order of clusters)

``` r
# order from heatmap
km_order <- row_order(Hmap)

# first cluster
cluster_1 <- as.data.frame(km_order[1]) %>% 
  mutate(k = "Cl 1") %>% 
  rename(n_row = 1)

# second cluster
cluster_2 <- as.data.frame(km_order[2]) %>% 
  mutate(k = "Cl 2") %>% 
  rename(n_row = 1)


# combine all cluster
clusters <- rbind(cluster_1, cluster_2) 


# join to standardized expression data 
clusters <- clusters %>% 
  arrange(-desc(n_row)) %>% 
  select(-n_row)

Cluster_data <- cbind(anova_sig %>% as.data.frame() %>%  rownames_to_column('Protein_ids'), clusters) 
```

### tidy anova data for supplementary tables

``` r
anova_sig_suppl <- read.delim("one_way_anova_perseus.txt", sep = "\t", header = T) %>% 
                   select(-starts_with("LFQ")) %>% 
                   filter(C..ANOVA.Significant == "+") %>% 
                   rename(Protein.IDs = T..Protein_IDs) %>% 
                   left_join(proteingroups %>% select(Protein.IDs, Gene.Symbol, Fasta.headers)) %>% 
                   left_join(Cluster_data %>% select(Protein_ids, k) %>% rename(Protein.IDs = Protein_ids, Cluster = k)) %>% 
                   arrange(-desc(Cluster)) %>% 
                   select(Cluster, Gene.Symbol, Protein.IDs, Fasta.headers, C..ANOVA.Significant, N...Log.ANOVA.p.value, T..Significant.pairs) %>% 
                   mutate(N...Log.ANOVA.p.value = 10^-N...Log.ANOVA.p.value) %>% 
                   rename(Genes                       = Gene.Symbol,
                          'Protein group'             = Protein.IDs,
                          'Fasta headers'             = Fasta.headers,
                          'One-way ANOVA significant (FDR<0.05)' = C..ANOVA.Significant,
                          'p-value'                   = N...Log.ANOVA.p.value,
                          'THSD pair'                 = T..Significant.pairs)
```


### get the order of samples on heatmap to make the same order in the profile plot

``` r
set.seed(1234)
animal_order <- column_order(Hmap) %>% 
                as.data.frame() %>% 
                rename(position=1) %>% 
                mutate(order = seq_along(1:ncol(anova_sig))) %>% 
                left_join(colnames(anova_sig) %>% 
                as.data.frame() %>% 
                rename(Animal = 1) %>% 
                mutate(position = seq_along(1:ncol(anova_sig)))) %>% 
                select(Animal) 

# facet labels
f_labels <- data.frame(Cluster = c("Cl 1", "Cl 2"), label = c("Cl 1", "Cl 2"))


# convert data to long table format
Cluster_data_long <- Cluster_data %>% 
  pivot_longer(names_to = "Bioreplicate", values_to = "Intensity", -c(Protein_ids, k))

# assign grouping
Cluster_data_long  <- Cluster_data_long %>% 
  left_join(conditions %>% 
              mutate(Bioreplicate = str_replace(Bioreplicate, "LFQ.intensity.", "")))

# make sure the order of facets corresponds order on heatmap
Cluster_data_long$k = factor(Cluster_data_long$k, levels= unique(anova_sig_suppl$Cluster)) 


# plot 
cluster_plot <- ggplot(Cluster_data_long, aes(x=Bioreplicate, y=as.numeric(Intensity))) +
  geom_line(aes(group=Protein_ids), size=0.4, alpha=0.5, color="grey")+
  geom_hline(yintercept=0, size= 0.5, linetype="dashed", color="firebrick3")+
  stat_summary(aes(group=1), fun=mean, color= "black", geom="line", size=0.6) +
  stat_summary(aes(group=1, color = Group),  shape = 15, fun=mean, geom="point", size=1.5)+
  scale_color_manual(values=c("Δ52" = "#708090", 'WT' = "navy", "Δ51-52" = "#D55E00")) +
  facet_grid(vars(k), scales = "free")+
  ylab("Intensity (z-scored)")      +
  scale_x_discrete(limits = animal_order$Animal) +
  xlab("") + 
  theme_bw() + 
  theme(plot.margin = margin(1,1,1,-0.5, "mm"))+
  theme(panel.border = element_rect(size = 1, colour = "black"),
        axis.title = element_text(size = 9, colour="black"), 
        axis.text.x = element_text(size= 9, colour="black", vjust = 0.5, angle = 90), 
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(size = 9, colour="black"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  geom_text(x = Inf, y = -Inf, aes(label = label),  size = 3, hjust = 1.05, vjust = -0.7, data = f_labels, inherit.aes = F)+
  theme(strip.background = element_blank(), strip.text = element_blank())
```

## over-representation analysis using webgestalt

### prepare the data for each cluster

``` r
# joined dataset
ora <- Cluster_data %>% 
  select(k, Protein_ids) %>% 
  rename(Protein.IDs = Protein_ids) %>% 
  left_join(proteingroups %>%  select(Gene.Symbol, Protein.IDs)) %>% 
  select(-Protein.IDs)
```


``` r
ora_for_cluster_function <- function(data, cluster_name) {
  
  # subset for specific cluster
  data_ora <- ora %>% 
  filter(k == cluster_name) %>% 
  select(-k) %>% 
  rename(Gene = 1) %>% 
  as.list()
  
  # perform ora
  set.seed(1234)
  outputDirectory <- getwd()
  enrichResult    <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
  enrichDatabase=c("pathway_KEGG", "geneontology_Biological_Process_noRedundant"), interestGene=data_ora[[1]],
  interestGeneType="genesymbol", referenceSet = "genome",
  referenceGeneType="genesymbol", isOutput=TRUE,
  outputDirectory=outputDirectory, projectName=paste0(cluster_name), sigMethod = "fdr", fdrThr = 0.25, fdrMethod = "BH", minNum = 5, maxNum = 200)
  
  # read ora output
  ora_results <- read.delim(paste0("Project_", cluster_name, "/", "enrichment_results_", cluster_name, ".txt"))
  
  # tidy ora output
  ora_results <- ora_results %>% 
    mutate(Cluster = cluster_name)
  
  return(ora_results)
}

# apply the function
cluster_1 <- ora_for_cluster_function(ora, cluster_name = "Cl 1")

cluster_2 <- ora_for_cluster_function(ora, cluster_name = "Cl 2")

# combine ora data
ora_data <- cluster_1 %>% 
  bind_rows(cluster_2) %>% 
  select(Cluster, geneSet, description, size, overlap, enrichmentRatio, pValue, FDR, userId) %>% 
  rename('id'                = geneSet, 
         'Biological process'   = description, 
         'Gene set size'        = size,
         'Gene count'           = overlap,
         'p-value'              = pValue,
         'Genes mapped'         = userId,
         'Fold enrichment'      = enrichmentRatio) %>% 
  mutate(`Biological process` = str_to_sentence(`Biological process`))

# save results
write.table(ora_data, "ora_data_each_cluster.txt", sep = "\t", row.names = F, quote = F)
rm(ora, km_order, animal_order, cluster_1, cluster_2)
```

### choose the processes that will be displayed on the plot

``` r
# replace zero fdr with lowest reported fdr (if any)
ora_data$FDR[ora_data$FDR == 0] <- min(ora_data$FDR[ora_data$FDR > 0])

### choose the processes that will be displayed on the plot
ora_data_plot <- ora_data 
```

### plot ora

``` r
# order rows based on enrichment score
ora_data_plot$`Biological process` <- factor(ora_data_plot$`Biological process`, levels = ora_data_plot$`Biological process`
                             [order(ora_data_plot$`Fold enrichment`)])

# plot
plot_dot <- ggplot(ora_data_plot, aes(x = `Fold enrichment`, y= `Biological process`, 
                                      fill = -log10(FDR), size = `Gene count`)) +
            geom_point(shape = 21)+
            theme_bw() +
            theme(panel.border = element_rect(size = 1, colour = "black"),
                            axis.text.x = element_text(angle = 90, colour = "black", vjust = 0.5, size = 8.5),
                            axis.title = element_text(size  = 9),
                            axis.text.y = element_text(size = 9, colour = "black"))+
            scale_x_continuous(breaks = c(20, 65,110))+
                      xlab("Enrichment") +
                      ylab("")+
                      theme(plot.title = element_text(size = 9, hjust=0.5,
                                                      face = "bold")) +
                      theme(plot.margin = margin(1,1,1,-1, "mm"))+
                      
                      theme(panel.grid.major = element_blank(),  
                            axis.ticks = element_line(colour = "black"),
                            panel.grid.minor = element_blank())+
                      scale_size_continuous(name = "Gene count", range = c(2,4), breaks = c(2,3,4))+
                      scale_fill_gradient(name = "-log10(FDR)", low = "#bf9deb", high = "#200b3b", breaks = c(1,2,3))+
                      theme(legend.position = "none", legend.key.size = unit(0.427, 'cm'), 
                            legend.box.spacing = unit(0.5, 'mm'), 
                            legend.title = element_text(size = 6), legend.text = element_text(size = 6))+
             facet_grid(vars(Cluster), scales = "free")+
             geom_text(x = Inf, y = -Inf, aes(label = label),  size = 3, hjust = 1.05, vjust = -0.7, data = f_labels, inherit.aes = F)+
             theme(strip.background = element_blank(), strip.text = element_blank())
```

### get legend from ORA (easier to arrange)

``` r
# fill legend
plot_dot_legend_fill <- plot_dot + 
             theme(legend.position = "bottom", legend.key.size = unit(0.427, 'cm'), 
                            legend.box.spacing = unit(0.5, 'mm'),
                            legend.title = element_text(size = 8.5), legend.text = element_text(size = 8.5))+
                            guides(size = "none")
# size legend
plot_dot_legend_size <- plot_dot + 
             theme(legend.position = "bottom", legend.key.size = unit(0.427, 'cm'), 
                            legend.box.spacing = unit(0.5, 'mm'),
                            legend.title = element_text(size = 8.5), legend.text = element_text(size = 8.5))+
                            guides(fill = "none")
# draw legends only (fill)
plot_dot_legend_fill <- get_legend(plot_dot_legend_fill) 
grid.newpage()                                     
grid.draw(plot_dot_legend_fill) 

# draw legends only (size)
plot_dot_legend_size <- get_legend(plot_dot_legend_size ) 
grid.newpage()                                     
grid.draw(plot_dot_legend_size ) 

# save legend
ggarrange(plot_dot_legend_fill, plot_dot_legend_size, nrow = 2)

ggsave("ora_legend.svg",  width =2.8, height = 0.9)
rm(plot_dot_legend_fill, plot_dot_legend_size)
```

## principal component analysis

### calculate principal components

``` r
pca_fn <- function(data, conditions_file) {
  
          # caclulate principal components of transposed dataframe
          PCA <- prcomp(t(data))
          
          # export results for ggplot
          # according to https://www.youtube.com/watch?v=0Jp4gsfOLMs (StatQuest: PCA in R)
          pca.data <- data.frame(Bioreplicate=rownames(PCA$x),
          X=PCA$x[,1],
          Y=PCA$x[,2])
          pca.var <- PCA$sdev^2
          pca.var.per <- round(pca.var/sum(pca.var)*100,1)
          # assign conditions
          pca.data <- pca.data %>% 
          left_join(conditions_file)
         
          # save the data in a list
          pca_data      <- list()
          pca_data[[1]] <- pca.data
          pca_data[[2]] <- pca.var.per[1:2]
         
          return(pca_data)
}
pca <- pca_fn(PG_imputed, conditions)
```

### PCA plot

``` r
# plot                     
pca_plot <- ggplot(data=pca[[1]], aes(x=X, y=Y, fill= Group))+
  geom_hline(yintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_point(size = 2.3, shape = 22) +
scale_fill_manual(values=c('Δ51-52' = "#D55E00", 'WT' = "navy", 'Δ52' = "#708090")) +
xlab(paste("PC 1 - ", pca[[2]][1], "%", sep=""))+
ylab(paste("PC 2 - ", pca[[2]][2], "%", sep=""))+
annotate("text", x = -14, y = 3, size = 3, label = "WT")+
annotate("text", x = 4, y = -12, size = 3, label = "Δ52")+
annotate("text", x = 0, y = 15, size = 3, label   = "Δ51-52")+
theme_bw() + 
  theme(panel.border = element_rect(size = 1, color = "black"),
        axis.title = element_text(size =  9, colour="black"),
        axis.ticks = element_line(colour = "black"),
        axis.text.x = element_text(size=  9, colour="black"), 
        axis.text.y = element_text(size = 9, colour="black"),
panel.grid.major = element_blank(), 
panel.grid.minor = element_blank())+
                      theme(legend.position = "top", 
                            legend.box.spacing = unit(0.5, 'mm'), 
                            legend.title = element_blank(), 
                            legend.text = element_text(size = 8.5),
                            legend.spacing.x  = unit(0.1, 'mm'),
                                 legend.margin=margin(0,0,0,0),
                            legend.box.margin=margin(-2.38,-2,-2.5,-2))
```

## volcano plots

### prepare data sets

``` r
volcano_fn <- function(data, xy, fc_relative_to) {
  
  # make an empty list
  data_volcano <- list()
  
  # tidy data for volcano
  data_tidy <-  data %>% 
    select(1:4) %>% 
    rename(Protein.IDs = Protein_IDs) %>% 
    left_join(proteingroups %>% select(Protein.IDs, Fasta.headers, Gene.Symbol)) %>% 
    select(Gene.Symbol, Protein.IDs, Fasta.headers, Significant, X.Log.P.value., Difference) %>% 
    mutate(X.Log.P.value. = 10^-X.Log.P.value.) %>% 
    rename(Gene            = Gene.Symbol,
           'Protein group' = Protein.IDs, 
           'Fasta headers' = Fasta.headers,
           'p-value'       = X.Log.P.value.,
           'log2 fold change'          = Difference) %>% 
    mutate(Diff_abundant = case_when(Significant == "+" & `log2 fold change` > 0 ~ paste0("increased_in_", fc_relative_to),
                                     Significant == "+" & `log2 fold change` < 0 ~ paste0("decreased_in_", fc_relative_to),
                                                              TRUE ~ "n.s.")) %>% 
    mutate_at(vars(contains("Significant")), ~replace(., is.na(.), "")) %>% 
    mutate(Significant = case_when(Significant == ""  ~ "n.s.",
                                   TRUE ~ Significant))
  
  
    
  # tidy s0 curve
   maxfc <- ceiling(max(data_tidy$`log2 fold change`))
   minfc <- floor(min(data_tidy$`log2 fold change`))
  
  if(maxfc > abs(minfc)) {
   
    z <- maxfc 
 
     } 
  
    else {z <- minfc}
  # filter curve
  Curve  <- xy %>% 
  filter(x > 0 & x < abs(z) | x < 0 & x > -abs(z)) 
  
  # return data
  data_volcano[[1]] <- data_tidy
  data_volcano[[2]] <- Curve
  
  return(data_volcano)
}
```

### volcano DMD vs WT

``` r
# apply the function
volcano_DMD <- volcano_fn(data = DMD_WT_perseus, xy = xy_DMD_WT, fc_relative_to = "DMD")

# plot volcano
volcano_plot_DMD <- ggplot(volcano_DMD[[1]] %>%                       
                       arrange(desc(`Diff_abundant`)),  mapping = aes(x = `log2 fold change`, y = -log10(`p-value`),
                                                           fill= Diff_abundant, 
                                                            alpha = Diff_abundant, 
                                                            shape = Diff_abundant)) +
                    geom_point(stroke = 0.2, size = 2.2) +
                    scale_fill_manual(values=c("n.s." = "#4a4949", "increased_in_DMD" = "firebrick3", "decreased_in_DMD" = "#56B4E9"))  +  
                    scale_shape_manual(values=c("n.s." = 16,       "increased_in_DMD" =21, "decreased_in_DMD" =21))+
                    scale_alpha_manual(values= c("n.s." = 0.3,     "increased_in_DMD" = 1, "decreased_in_DMD" = 1)) +
  geom_path(data = volcano_DMD[[2]], mapping = aes(x ,y), size =0.5 ,inherit.aes = FALSE, na.rm = TRUE, color = "black") +
  theme(panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.background = element_blank(), 
       axis.line = element_line(color = "black"),
  legend.position = "NONE") +
  theme(axis.title = element_text(size = 8.5), 
        axis.text.x = element_text(size= 8.5, colour = "black", vjust = -0.1), 
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(size = 8.5, colour = "black"))+
  xlab("log2 fold change (Δ52/WT)")+
  ylab("-log10 p-value") +
  scale_y_continuous(breaks = c(0,2,4,6,8),   limits = c(0,   9))+
  scale_x_continuous(breaks = c(-8,-4,0,4,8), limits = c(-8.6, 8.6))+
  annotate("text", x =  7, y = 9, size  = 3.3, label = paste("↑",  round(((sum(str_count(volcano_DMD[[1]], "increased_in_DMD")))/nrow(volcano_DMD[[1]]))*100, digits = 1), "%"))+
  annotate("text", x = -6, y = 9, size = 3.3, label = paste("↓", round(((sum(str_count(volcano_DMD[[1]], "decreased_in_DMD")))/nrow(volcano_DMD[[1]]))*100, digits = 1), "%"))+
  ggtitle("Δ52") +                
  theme(plot.title = element_text(hjust = 0.5, size = 8.5))
```


### volcano BMD vs WT

``` r
# apply the function
volcano_BMD <- volcano_fn(data = BMD_WT_perseus, xy = xy_BMD_WT, fc_relative_to = "BMD")

# plot volcano
volcano_plot_BMD <- ggplot(volcano_BMD[[1]] %>%                       
                       arrange(desc(`Diff_abundant`)),  mapping = aes(x = `log2 fold change`, y = -log10(`p-value`),
                                                            fill= `Diff_abundant`, 
                                                            alpha = `Diff_abundant`, 
                                                            shape = `Diff_abundant`)) +
                    geom_point(stroke = 0.2, size = 2.2) +
                    scale_fill_manual(values=c("n.s." = "#4a4949", "increased_in_BMD" = "firebrick3", "decreased_in_BMD" = "#56B4E9"))  +  
                    scale_shape_manual(values=c("n.s." = 16, "increased_in_BMD" =21, "decreased_in_BMD" =21))+
                    scale_alpha_manual(values= c("n.s." = 0.3, "increased_in_BMD" = 1, "decreased_in_BMD" = 1)) +
  geom_path(data = volcano_BMD[[2]], mapping = aes(x ,y), size =0.5 ,inherit.aes = FALSE, na.rm = TRUE, color = "black") +
  theme(panel.grid.major = element_blank(), 
       panel.grid.minor = element_blank(), 
       panel.background = element_blank(), 
       axis.line = element_line(color = "black"),
  legend.position = "NONE") +
  theme(axis.title = element_text(size = 8.5), 
        axis.text.x = element_text(size= 8.5, colour = "black", vjust = -0.1), 
        axis.ticks = element_line(colour = "black"),
        axis.text.y = element_text(size = 8.5, colour = "black"))+
  xlab("log2 fold change (Δ51-52/WT)")+
  ylab("-log10 p-value") +
  scale_y_continuous(breaks = c(0,2,4,6), limits = c(0,6))+
  scale_x_continuous(breaks = c(-6,-3,0,3,6))+
  annotate("text", x = 4, y = 6, size = 3, label = paste("↑",  round(((sum(str_count(volcano_BMD[[1]], "increased_in_BMD")))/nrow(volcano_BMD[[1]]))*100, digits = 1), "%"))+
  annotate("text", x = -5, y = 6, size = 3, label = paste("↓", round(((sum(str_count(volcano_BMD[[1]], "decreased_in_BMD")))/nrow(volcano_BMD[[1]]))*100, digits = 1), "%"))+
  ggtitle("Δ51-52") +                
  theme(plot.title = element_text(hjust = 0.5, size = 8.5))
```


## save data

### combine plots (main figure)

``` r
p1 <- ggarrange(pca_plot, volcano_plot_DMD, volcano_plot_BMD, labels = c("c", "d", "e"), font.label = list(size = 14, face = "bold"),
          ncol = 3, nrow = 1, widths = c(7.1/3,7.1/3,7.1/3), heights = c(7.1/3,7.1/3,7.1/3))
p2 <- rectGrob(width = 1, height = 1)
p3 <- ggarrange(cluster_plot, p2, widths = c(1.8,1.8), heights = c(h1-1, 1), ncol = 1, labels = c("g"), font.label = list(size = 14, face = "bold"))
p4 <- ggarrange(heatmap_plot, p3, labels = c("f"), font.label = list(size = 14), widths = c(w1, 2))
p5 <- ggarrange(plot_dot, p2, widths = c(7.1-(w1+2),7.1-(w1+2)), heights = c(h1-1, 1), ncol = 1, labels = c("h"), font.label = list(size = 14, face = "bold"))
p6 <- ggarrange(p4, p5, heights = c(h1, h1), widths = c(w1+2, 7.1-(w1+2)))
p7 <- ggarrange(p1, p6, ncol = 1, widths = c(w1+2+7.1-(w1+2), w1+2+7.1-(w1+2)), heights = c(2.1, h1))
ggsave("BMD_myoc.svg", width = w1+2+7.1-(w1+2), height = h1+2.1)
rm(p1,p2,p3,p4,p5,p6,p7)
```

### tidy proteingroups file for supplementary tables

``` r
# load protein groups file again without formating column names
rm(proteingroups)
proteingroups <- read.delim("proteinGroups_myocardium.txt", sep = "\t", header = T, check.names = F)
# tidy for supplementary table
PG_suppl <- proteingroups %>% 
  select(`Gene Symbol`,
         1:2,
         `Fasta headers`,
         `Unique peptides`,
         `Razor + unique peptides`,
         contains("Intensity"), 
         contains("MS/MS count"), 
         `Only identified by site`, 
         Reverse, 
         `Potential contaminant`, 
         contains("Sequence coverage [%]"),
         `Mol. weight [kDa]`
) %>% 
  arrange(desc(Intensity))
```

### save data in a supplementary excel files

``` r
# protein groups
write.xlsx(PG_suppl, file = "proteomics_suppl.xlsx", sheetName = "Suppl table 1A", 
  col.names = TRUE, row.names = FALSE, append = T)
# anova results 
write.xlsx(anova_sig_suppl, file = "proteomics_suppl.xlsx", sheetName = "Suppl table 1B", 
  col.names = TRUE, row.names = FALSE, append = T)
# ora data
write.xlsx(ora_data, file = "proteomics_suppl.xlsx", sheetName = "Suppl table 1C", 
  col.names = TRUE, row.names = FALSE, append = T)
```
