library(phyloseq)
library(ggplot2)
library(qiime2R)
library(data.table)
library(DESeq2)
library(dplyr)
library(edgeR)

myTheme <- theme(panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 text = element_text(color = 'black', size = 12),
                 axis.text = element_text(color = 'black'))

# Note that for proper loading into R, you need to delete the # character before 'SampleID'
# In the metadata txt file

CPseq <- qza_to_phyloseq(features = "protein_qiime2_qza_files/CP_table.qza",
                         taxonomy = "protein_qiime2_qza_files/CP_taxonomy.qza",
                         tree = "protein_qiime2_qza_files/CP_rooted_tree.qza",
                         metadata = "protein_qiime2_featuretables/CP_mapping_file.txt")

#CPseq <- subset_samples(CPseq, Timepoint != "D-1")
CPseq <- subset_samples(CPseq, Timepoint %in% c("D0", "D4", "D8"))
CPseq <- subset_samples(CPseq, Treatment_Group2 %in% c("pct00", "pct06", "pct18"))
CPseq@sam_data$Treatment_Group2 <- factor(CPseq@sam_data$Treatment_Group2,
                                          labels = c("0%", "6%", "18%"))

GPseq <- qza_to_phyloseq(features = "protein_qiime2_qza_files/GP_table.qza",
                         taxonomy = "protein_qiime2_qza_files/GP_taxonomy.qza",
                         tree = "protein_qiime2_qza_files/GP_rooted_tree.qza",
                         metadata = "protein_qiime2_featuretables/GP_mapping.txt")

GPseq <- subset_samples(GPseq, Timepoint %in% c("D07", "D14", "D19"))
GPseq@sam_data$Treatment_Group2 <- factor(GPseq@sam_data$Treatment_Group2,
                                          labels = c("0%", "6%", "18%"))

PTseq <- qza_to_phyloseq(features = "protein_qiime2_qza_files/PT_table.qza",
                         taxonomy = "protein_qiime2_qza_files/PT_taxonomy.qza",
                         tree = "protein_qiime2_qza_files/PT_rooted_tree.qza",
                         metadata = "protein_qiime2_featuretables/PT_N_L_mapping.txt")

PTseq <- subset_samples(PTseq, Timepoint %in% c("D-2", "D0", "D1", "D2", "D3", "D7"))

# check the ranks and samples to ensure proper loading
rank_names(CPseq)
rank_names(PTseq)
rank_names(GPseq)

sample_names(CPseq)
sample_names(PTseq)
sample_names(GPseq)

###################################
# CP Analysis
###################################
CPseq <- subset_taxa(CPseq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

CPprev <- apply(X = otu_table(CPseq),
                MARGIN = ifelse(taxa_are_rows(CPseq), yes = 1, no = 2),
                FUN = function(x) {sum(x > 0)})

CPprev <- data.frame(Prevalence = CPprev,
                     TotalAbundance = taxa_sums(CPseq),
                     tax_table(CPseq))

ggplot(CPprev, aes(x = TotalAbundance, y = Prevalence / nsamples(CPseq), color = Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +labs(x = "Total Abundance", y = "Prevalence (fraction of samples)") +
  facet_wrap(~Phylum, 3) + theme(legend.position = "none") + myTheme

plyr::ddply(CPprev, "Phylum", function(df1)
  {cbind(mean(df1$Prevalence), sum(df1$Prevalence))})

CPseqFilter <- genefilter_sample(CPseq, filterfun_sample(function(x) x > 5), 
                                 A = 0.25*nsamples(CPseq))
CPseq <- prune_taxa(CPseqFilter, CPseq)

CPseq <- transform_sample_counts(CPseq, function(x) 100 * x/sum(x))
CPphylumSum <- tapply(taxa_sums(CPseq), tax_table(CPseq)[, "Phylum"], sum, na.rm = T)
CPtop5 <- names(sort(CPphylumSum, T))[1:5]
CPseq <- prune_taxa((tax_table(CPseq)[, "Phylum"] %in% CPtop5), CPseq)

CPord <- ordinate(CPseq, "NMDS", "bray")
#plot_ordination(CPseq, CPord, type = "taxa", color = "Phylum") +
#  facet_wrap(~Phylum)

ordPlot <- plot_ordination(CPseq, CPord, type = "samples", color = "Treatment_Group2") +
  facet_wrap(~Timepoint) + myTheme + #stat_ellipse(level = 0.8) +
  labs(color = "Percent protein") +
  scale_color_manual(values = c("#d73027", "#fee090", "#4575b4", "gray", "black")) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1),
        legend.key=element_blank())

library(vegan)
cpMeta <- as(sample_data(CPseq), "data.frame")
adonis(phyloseq::distance(subset_samples(CPseq, Timepoint == "D8"), 
                          method = "bray") ~ Treatment_Group2,
       data = cpMeta[cpMeta$Timepoint == "D8",])

ordPlotDat <- ordPlot$data
ggplot(ordPlotDat, aes(x = Treatment_Group2, y = NMDS1)) +
  geom_boxplot() + facet_wrap(~Timepoint) +
  labs(x = "Percent protein") + ggpubr::stat_compare_means(label.y.npc = 0.93, label.x.npc = 0.05) +
  myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1),
        legend.background=element_blank())

plot_richness(CPseq, measures = c("Shannon"),
             x = "Treatment_Group2") +
  facet_wrap(~Timepoint) + geom_boxplot() + myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = "Percent protein", y = "Shannon alpha diversity") +
  ggpubr::stat_compare_means(label.y.npc = 0.99, label.x.npc = 0.05)

plot_richness(CPseq, measures = c("Shannon"),
              x = "TimeNum", color = "Treatment_Group2") +
  stat_summary(aes(y = value), fun.y=mean, geom="line", size = 1) +
  scale_color_manual(values = c("#d73027", "#fee090", "#4575b4", "gray", "black")) +
  geom_point(size = 3) + myTheme + theme(axis.text.x = element_text(angle = 0)) +
  labs(x = "Days on diet", y = "Shannon alpha diversity", color = "Percent protein")

PCdistBC <- distance(CPseq, method = "bray")
PCdistUF <- distance(CPseq, method = "wUniFrac")

PCordBC <- ordinate(CPseq, method = "PCoA", distance = PCdistBC)
PCordUF <- ordinate(CPseq, method = "PCoA", distance = PCdistUF)

plot_scree(PCordBC)
plot_scree(PCordUF)

plot_ordination(CPseq, PCordBC, color = "Treatment_Group2") +
  geom_point() + facet_wrap(~Timepoint)

CPseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, Phylum) %>%
  summarise(Proportion = sum(Abundance, na.rm = T) / length(unique(MouseID))) %>%
  ggplot(aes(x = Treatment_Group2, y = Proportion, fill = Phylum)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  facet_wrap(~Timepoint) + 
  labs(x = "Percent protein", y = "Relative abundance") +
  myTheme + theme(panel.border = element_rect(color = "black", fill = NA),
                 strip.background = element_rect(color = "black", size = 1))

phylumFrame <- CPseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, MouseID, Phylum) %>%
  summarise(Proportion = sum(Abundance, na.rm = T))

genusFrame <- CPseq %>%
  psmelt %>%
  group_by(Treatment_Group2, Timepoint, MouseID, Genus) %>%
  summarise(Proportion = sum(Abundance, na.rm = T))

ggplot(phylumFrame, aes(x = Treatment_Group2, y = Proportion)) +
  facet_grid(Phylum ~ Timepoint, scales = "free") + geom_point() +
  ggpubr::stat_compare_means(label.y.npc = 0.9, label.x.npc = 0.1) + myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, aes(group = 1)) +
  labs(x = "Percent protein", y = "Relative abundance (%)")

# Oscillospira, Coprococcus
oscilloFrame <- subset(genusFrame, (!is.na(Treatment_Group2)) &
                         Genus == "Turicibacter")

ggplot(oscilloFrame, aes(x = Timepoint, y = Proportion)) +
  geom_point(shape = 21, position = position_jitterdodge(), 
             aes(fill = Treatment_Group2), color = "black") + 
  geom_boxplot(color = "black", aes(fill = Treatment_Group2)) + myTheme +
  scale_color_manual(values = c("#d73027", "#fee090", "#4575b4", "gray", "black")) +
  scale_fill_manual(values = c("#d73027", "#fee090", "#4575b4", "gray", "black")) +
  labs(y = "Turicibacter abundance (%)", color = "Percent protein", fill = "Percent protein") +
  theme(legend.key=element_blank())

## Differential Expression
## Note: Rerun phyloseq object generation and but do not relabel treatment factor
CPmat <- as(otu_table(CPseq), "matrix")
CPmat <- CPmat + 1

CPtaxonomy <- tax_table(CPseq, errorIfNULL = F)
if( !is.null(CPtaxonomy)){
  CPtaxonomy <- data.frame(as(CPtaxonomy, "matrix"))
}

CPdge <- DGEList(counts = CPmat, genes = CPtaxonomy, remove.zeros = T)
CPdge <- calcNormFactors(CPdge, method = "RLE")

plotMDS(CPdge)

CPdays <- CPseq@sam_data$Timepoint
CPprot <- CPseq@sam_data$Treatment_Group2

# Approach 1 - combine time and tx factors to explicitly specify contrasts
CPmm <- model.matrix(~ 0 + interaction(Treatment_Group2, Timepoint),
                     data = data.frame(CPseq@sam_data))
colnames(CPmm) <- stringr::str_sub(colnames(CPmm), -8, -1)

CPv <- voom(CPdge, CPmm, plot = T)
CPfit <- lmFit(CPv, CPmm)
colnames(coef(CPfit))
CPcontrasts <- makeContrasts(
  d0_18v0 = pct00.D0-pct18.D0,
  d0_18v6 = pct06.D0-pct18.D0,
  #d1_18v0 = pct00.D1-pct18.D1,
  #d1_18v6 = pct06.D1-pct18.D1,
  #d2_18v0 = pct00.D2-pct18.D2,
  #d2_18v6 = pct06.D2-pct18.D2,
  #d3_18v0 = pct00.D3-pct18.D3,
  #d3_18v6 = pct06.D3-pct18.D3,
  d4_18v0 = pct00.D4-pct18.D4,
  d4_18v6 = pct06.D4-pct18.D4,
  d8_18v0 = pct00.D8-pct18.D8,
  d8_18v6 = pct06.D8-pct18.D8,
  levels = CPmm
)

CPtmp <- contrasts.fit(CPfit, CPcontrasts)
CPtmp <- eBayes(CPtmp)
colnames(coef(CPtmp))

CPtt <- topTable(CPtmp, coef = 6, sort.by = "P", n = Inf)
CPtt$Taxa <- rownames(CPtt)
length(which(CPtt$adj.P.Val < 0.05))

CPsig <- cbind(as(CPtt, "data.frame"),
               as(tax_table(CPseq)[rownames(CPtt), ], "matrix"))

CPsig <- subset(CPsig, !is.na(Genus))

x = tapply(CPsig$logFC, CPsig$Phylum, function(x) max(x))
x = sort(x, T)
CPsig$Phylum <- factor(as.character(CPsig$Phylum), levels = names(x))

x = tapply(CPsig$logFC, CPsig$Genus, function(x) max(x))
x = sort(x, T)
CPsig$Genus <- factor(as.character(CPsig$Genus), levels = names(x))

CPsig$sig <- F
CPsig[CPsig$adj.P.Val < 0.05,]$sig <- T

CPsig$Phylum <- factor(CPsig$Phylum,
                       levels = c("Actinobacteria", "Bacteroidetes",
                                  "Firmicutes", "Proteobacteria",
                                  "Verrucomicrobia"))
  
ggplot(CPsig, aes(x = Genus, y = logFC)) +
  geom_point(shape = 21, aes(fill = Phylum, color = sig), size = 2.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = c("white", "black")) +
  labs( y = "Log2 fold change", color = "Significant", title = "6% v 18% D8") +
  myTheme + theme(legend.key=element_blank())

###################################
# GP Analysis
###################################
GPseq <- subset_taxa(GPseq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

GPprev <- apply(X = otu_table(GPseq),
                MARGIN = ifelse(taxa_are_rows(GPseq), yes = 1, no = 2),
                FUN = function(x) {sum(x > 0)})

GPprev <- data.frame(Prevalence = GPprev,
                     TotalAbundance = taxa_sums(GPseq),
                     tax_table(GPseq))

ggplot(GPprev, aes(x = TotalAbundance, y = Prevalence / nsamples(GPseq), color = Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +labs(x = "Total Abundance", y = "Prevalence (fraction of samples)") +
  facet_wrap(~Phylum, 3) + theme(legend.position = "none") + myTheme

plyr::ddply(GPprev, "Phylum", function(df1)
{cbind(mean(df1$Prevalence), sum(df1$Prevalence))})

GPseqFilter <- genefilter_sample(GPseq, filterfun_sample(function(x) x > 5), 
                                 A = 0.25*nsamples(GPseq))
GPseq <- prune_taxa(GPseqFilter, GPseq)

GPseq <- transform_sample_counts(GPseq, function(x) 100 * x/sum(x))
GPphylumSum <- tapply(taxa_sums(GPseq), tax_table(GPseq)[, "Phylum"], sum, na.rm = T)
GPtop5 <- names(sort(GPphylumSum, T))[1:5]
GPseq <- prune_taxa((tax_table(GPseq)[, "Phylum"] %in% GPtop5), GPseq)

GPord <- ordinate(GPseq, "NMDS", "bray")
plot_ordination(GPseq, GPord, type = "taxa", color = "Phylum") +
  facet_wrap(~Phylum)

GPordPlot <- plot_ordination(GPseq, GPord, type = "samples", color = "Treatment_Group2") +
  facet_wrap(~Timepoint) + myTheme + #stat_ellipse(level = 0.8) +
  labs(color = "% Protein") +
  scale_color_manual(values = c("#d73027", "#fee090", "#4575b4")) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1),
        legend.key=element_blank())

gpMeta <- as(sample_data(GPseq), "data.frame")
adonis(phyloseq::distance(subset_samples(GPseq, Timepoint == "D19"), 
                          method = "bray") ~ Treatment_Group2,
       data = gpMeta[gpMeta$Timepoint == "D19",])

GPordPlotDat <- GPordPlot$data
ggplot(GPordPlotDat, aes(x = Treatment_Group2, y = NMDS1)) +
  geom_boxplot() + facet_wrap(~Timepoint) +
  labs(x = "Percent protein") + ggpubr::stat_compare_means(label.y.npc = 0.93) +
  myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1))

plot_richness(GPseq, measures = c("Shannon"),
              x = "Treatment_Group2") +
  facet_wrap(~Timepoint) + geom_boxplot() + myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = "Percent protein", y = "Shannon alpha diversity") +
  ggpubr::stat_compare_means(label.y.npc = 0.99)

plot_richness(GPseq, measures = c("Shannon"),
              x = "TimeNum", color = "Treatment_Group2") +
  stat_summary(aes(y = value), fun.y=mean, geom="line", size = 1) +
  scale_color_manual(values = c("#d73027", "#fee090", "#4575b4")) +
  geom_point(size = 3) + myTheme + theme(axis.text.x = element_text(angle = 0)) +
  labs(x = "Days on diet", y = "Shannon alpha diversity", color = "Percent protein")

PCdistBC <- distance(GPseq, method = "bray")
PCdistUF <- distance(GPseq, method = "wUniFrac")

PCordBC <- ordinate(GPseq, method = "PCoA", distance = PCdistBC)
PCordUF <- ordinate(GPseq, method = "PCoA", distance = PCdistUF)

plot_scree(PCordBC)
plot_scree(PCordUF)

plot_ordination(GPseq, PCordBC, color = "Treatment_Group2") +
  geom_point() + facet_wrap(~Timepoint)

GPseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, Phylum) %>%
  summarise(Proportion = sum(Abundance, na.rm = T) / length(unique(MouseID))) %>%
  ggplot(aes(x = Treatment_Group2, y = Proportion, fill = Phylum)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  facet_wrap(~Timepoint) + 
  labs(x = "Percent protein", y = "Relative abundance") +
  myTheme + theme(panel.border = element_rect(color = "black", fill = NA),
                  strip.background = element_rect(color = "black", size = 1))

phylumFrameGP <- GPseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, MouseID, Phylum) %>%
  summarise(Proportion = sum(Abundance, na.rm = T))

genusFrameGP <- GPseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, MouseID, Genus) %>%
  summarise(Proportion = sum(Abundance, na.rm = T))

ggplot(phylumFrameGP, aes(x = Treatment_Group2, y = Proportion)) +
  facet_grid(Phylum ~ Timepoint, scales = "free") + geom_point() +
  ggpubr::stat_compare_means(label.y.npc = 0.9, label.x.npc = 0.1) + myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, aes(group = 1)) +
  labs(x = "Percent protein", y = "Relative abundance (%)")

oscilloFrameGP <- subset(genusFrameGP, (!is.na(Genus)) &
                          (Genus == "Turicibacter"))

ggplot(oscilloFrameGP, aes(x = Timepoint, y = Proportion)) +
  geom_point(shape = 21, position = position_jitterdodge(), 
             aes(fill = Treatment_Group2), color = "black") + 
  geom_boxplot(color = "black", aes(fill = Treatment_Group2)) + myTheme +
  scale_color_manual(values = c("#d73027", "#fee090", "#4575b4")) +
  scale_fill_manual(values = c("#d73027", "#fee090", "#4575b4")) +
  labs(y = "Turicibacter abundance (%)", color = "Percent protein", fill = "Percent protein") +
  theme(legend.key=element_blank())

## Differential Expression
GPmat <- as(otu_table(GPseq), "matrix")
GPmat <- GPmat + 1

GPtaxonomy <- tax_table(GPseq, errorIfNULL = F)
if( !is.null(GPtaxonomy)){
  GPtaxonomy <- data.frame(as(GPtaxonomy, "matrix"))
}

GPdge <- DGEList(counts = GPmat, genes = GPtaxonomy, remove.zeros = T)
GPdge <- calcNormFactors(GPdge, method = "RLE")

plotMDS(GPdge)

GPdays <- GPseq@sam_data$Timepoint
GPprot <- GPseq@sam_data$Treatment_Group2

# Approach 1 - combine time and tx factors to explicitly specify contrasts
GPmm <- model.matrix(~ 0 + interaction(Treatment_Group2, Timepoint),
                     data = data.frame(GPseq@sam_data))
colnames(GPmm) <- stringr::str_sub(colnames(GPmm), -9, -1)
colnames(GPmm) <- gsub(")", "", colnames(GPmm))


GPv <- voom(GPdge, GPmm, plot = T)
GPfit <- lmFit(GPv, GPmm)
colnames(coef(GPfit))
GPcontrasts <- makeContrasts(
  #d0_18v0 = pct00.D00-pct18.D00,
  #d0_18v6 = pct06.D00-pct18.D00,
  d7_18v0 = pct00.D07-pct18.D07,
  d7_18v6 = pct06.D07-pct18.D07,
  d14_18v0 = pct00.D14-pct18.D14,
  d14_18v6 = pct06.D14-pct18.D14,
  d21_18v0 = pct00.D19-pct18.D19,
  d21_18v6 = pct06.D19-pct18.D19,
  #dbase_18v0 = pct00.IN-pct18.IN,
  #dbase_18v6 = pct06.IN-pct18.IN,
  levels = GPmm
)

GPtmp <- contrasts.fit(GPfit, GPcontrasts)
GPtmp <- eBayes(GPtmp)
colnames(coef(GPtmp))

GPtt <- topTable(GPtmp, coef = 6, sort.by = "P", n = Inf)
GPtt$Taxa <- rownames(GPtt)
length(which(GPtt$adj.P.Val < 0.05))

GPsig <- cbind(as(GPtt, "data.frame"),
               as(tax_table(GPseq)[rownames(GPtt), ], "matrix"))

GPsig <- subset(GPsig, !is.na(Genus))

x = tapply(GPsig$logFC, GPsig$Phylum, function(x) max(x))
x = sort(x, T)
GPsig$Phylum <- factor(as.character(GPsig$Phylum), levels = names(x))

x = tapply(GPsig$logFC, GPsig$Genus, function(x) max(x))
x = sort(x, T)
GPsig$Genus <- factor(as.character(GPsig$Genus), levels = names(x))

GPsig$sig <- F
GPsig[GPsig$adj.P.Val < 0.05,]$sig <- T

GPsig$Phylum <- factor(GPsig$Phylum,
                       levels = c("Actinobacteria", "Bacteroidetes",
                                  "Firmicutes", "Proteobacteria",
                                  "Verrucomicrobia"))

ggplot(GPsig, aes(x = Genus, y = logFC)) +
  geom_point(shape = 21, aes(fill = Phylum, color = sig), size = 2.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_color_manual(values = c("white", "black")) +
  labs( y = "Log2 fold change", color = "Significant", title = "6% v 18% D19") +
  myTheme + theme(legend.key=element_blank())

##########################
# Combined CP GP EdgeR
##########################

CPGPedge <- merge(CPtt, GPtt, by = "Taxa", suffixes = c(".CP", ".GP"))
CPGPedge <- CPGPedge[!is.na(CPGPedge$Genus.CP), ]

CPGPedge$logFC.CP
ggplot(CPGPedge, aes(x = logFC.CP, y = logFC.GP, color = Genus.CP)) +
  geom_point() + myTheme +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Donor log FC", y = "Recipient log FC", color = "Genus") +
  xlim(0, 1.5) + ylim(0, 0.6)

###################################
# PT Analysis
###################################
PTseq <- subset_taxa(PTseq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

PTprev <- apply(X = otu_table(PTseq),
                MARGIN = ifelse(taxa_are_rows(PTseq), yes = 1, no = 2),
                FUN = function(x) {sum(x > 0)})

PTprev <- data.frame(Prevalence = PTprev,
                     TotalAbundance = taxa_sums(PTseq),
                     tax_table(PTseq))

ggplot(PTprev, aes(x = TotalAbundance, y = Prevalence / nsamples(PTseq), color = Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +labs(x = "Total Abundance", y = "Prevalence (fraction of samples)") +
  facet_wrap(~Phylum, 3) + theme(legend.position = "none") + myTheme

plyr::ddply(PTprev, "Phylum", function(df1)
{cbind(mean(df1$Prevalence), sum(df1$Prevalence))})

PTseqFilter <- genefilter_sample(PTseq, filterfun_sample(function(x) x > 5), 
                                 A = 0.25*nsamples(PTseq))
PTseq <- prune_taxa(PTseqFilter, PTseq)

PTseq <- transform_sample_counts(PTseq, function(x) 100 * x/sum(x))
PTphylumSum <- tapply(taxa_sums(PTseq), tax_table(PTseq)[, "Phylum"], sum, na.rm = T)
PTtop5 <- names(sort(PTphylumSum, T))[1:5]
PTseq <- prune_taxa((tax_table(PTseq)[, "Phylum"] %in% PTtop5), PTseq)

PTord <- ordinate(PTseq, "NMDS", "bray")
plot_ordination(PTseq, PTord, type = "taxa", color = "Phylum") +
  facet_wrap(~Phylum)
plot_ordination(PTseq, PTord, type = "samples", color = "Treatment_Group2") +
  facet_wrap(~Timepoint) + myTheme + #stat_ellipse(level = 0.8) +
  labs(color = "% Protein") +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1))

plot_richness(PTseq, measures = c("Shannon"),
              x = "Treatment_Group2") +
  facet_wrap(~Timepoint) + geom_boxplot() + myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1),
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  labs(x = "Percent protein", y = "Shannon alpha diversity") +
  ggpubr::stat_compare_means()

plot_richness(PTseq, measures = c("Shannon"),
              x = "TimeNum", color = "Treatment_Group2") +
  stat_summary(aes(y = value), fun.y=mean, geom="line", size = 1) +
  geom_point(size = 3) + myTheme + theme(axis.text.x = element_text(angle = 0)) +
  labs(x = "Days on diet", y = "Shannon alpha diversity", color = "% protein")

PCdistBC <- distance(PTseq, method = "bray")
PCdistUF <- distance(PTseq, method = "wUniFrac")

PCordBC <- ordinate(PTseq, method = "PCoA", distance = PCdistBC)
PCordUF <- ordinate(PTseq, method = "PCoA", distance = PCdistUF)

plot_scree(PCordBC)
plot_scree(PCordUF)

plot_ordination(PTseq, PCordBC, color = "Treatment_Group2") +
  geom_point() + facet_wrap(~Timepoint)

PTseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, Phylum) %>%
  summarise(Proportion = sum(Abundance, na.rm = T) / length(unique(MouseID))) %>%
  ggplot(aes(x = Treatment_Group2, y = Proportion, fill = Phylum)) +
  geom_bar(stat = 'identity', position = 'stack', color = 'black') +
  facet_wrap(~Timepoint) + 
  labs(x = "Percent protein", y = "Relative abundance") +
  myTheme + theme(panel.border = element_rect(color = "black", fill = NA),
                  strip.background = element_rect(color = "black", size = 1))

phylumFrame <- PTseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, MouseID, Phylum) %>%
  summarise(Proportion = sum(Abundance, na.rm = T))

ggplot(phylumFrame, aes(x = Treatment_Group2, y = Proportion)) +
  facet_grid(Phylum ~ Timepoint, scales = "free") + geom_point() +
  ggpubr::stat_compare_means(label.y.npc = 0.9, label.x.npc = 0.1) + myTheme +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_rect(color = "black", size = 1)) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x, aes(group = 1)) +
  labs(x = "Percent protein", y = "Relative abundance (%)")


genusFramePT <- PTseq %>%
  psmelt() %>%
  group_by(Treatment_Group2, Timepoint, MouseID, Genus) %>%
  summarise(Proportion = sum(Abundance, na.rm = T))

oscilloFramePT <- subset(genusFramePT, (!is.na(Genus)) &
                           (Genus == "Oscillospira"))

oscilloFramePT$timeNum <- (factor(oscilloFramePT$Timepoint,
                                  labels = c("-2", "0", "1", "2", "3", "7")))

oscilloFramePT$timeNum <- as.numeric(as.character(oscilloFramePT$timeNum))

ggplot(oscilloFramePT, aes(x = timeNum, y = Proportion, color = Treatment_Group2)) +
  geom_point() + geom_smooth(se = F)

###################################
# PTCP Analysis
###################################
PTCPseqFilter <- genefilter_sample(PTCPseq, filterfun_sample(function(x) x > 5), 
                                   A = 0.5*nsamples(PTCPseq))
PTCPseq <- prune_taxa(PTCPseqFilter, PTCPseq)

PTCPseq <- transform_sample_counts(PTCPseq, function(x) 1E6 * x/sum(x))
PTCPphylumSum <- tapply(taxa_sums(PTCPseq), tax_table(PTCPseq)[, "Phylum"], sum, na.rm = T)
PTCPtop5 <- names(sort(PTCPphylumSum, T))[1:5]
PTCPseq <- prune_taxa((tax_table(PTCPseq)[, "Phylum"] %in% PTCPtop5), PTCPseq)

PTCPord <- ordinate(PTCPseq, "NMDS", "bray")
plot_ordination(PTCPseq, PTCPord, type = "taxa", color = "Phylum") +
  facet_wrap(~Phylum)
plot_ordination(PTCPseq, PTCPord, type = "samples", color = "Treatment_Group2") +
  facet_wrap(~Timepoint)

plot_richness(PTCPseq, measures = c("Shannon"),
              x = "Treatment_Group2", color = "Treatment_Group2") +
  facet_wrap(~Timepoint)
