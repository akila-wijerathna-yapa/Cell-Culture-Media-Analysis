#===============================================================================
#Combe 2025 for "Comparative Metabolomic Profiling of Hydrolysates from Yeast -
#   and Plant-Based Sources for Cell Culture Media Development" 
#-------------------------------------------------------------------------------
library(cowplot)
library(ggrepel)
library(tidyverse)
library(VennDiagram)
library(pheatmap)

#===============================================================================
#Plot colour options

Col_all = c("Pea" = "#9f79EE","Soy" =  "#EE6363", 
            "Wheat" = "#7EC0EE", "Cotton" = "#FFC125",
            "YE" = "#C3D7A4", "Hy-Yest-466" = "#4EEE94", 
            "Hy-Yest-412" = "#52854C",
            "Hy-Yest-503" = "#7CCD7C", "Hy-Yest-555" ="#8F9779")

Col_metab = c("#BFD3C1", "#F4A6A1", "#C9DAF8", "#FAD02E", "#D5A6BD",
              "#A3D9A1", "#F4B6C2", "#B0E0E6", "#F5C6CB", "#D6E6F2",
              "#F9E2AE", "#C2C2F0", "#E1BEE7", "#A7D4C6", "#F8C3A4")

Col2 = c("#9f79EE", "#EE6363", "#CD661D","#7EC0EE","#FFC125", "#52854C",
         "#1C86EE", "#7CCD7C", "#8B0A50")

col_cat = c( "Amino Acid" = "#7AC5CD","Carbohydrates" = "#FFABA9",
             "Lipids" = "#7CCD7E","Nucleosides" = "#AB82FF", 
             "Organic acids" = "#FFA54F", "Other" = "#EEAEEE", 
             "Peptides" = "#FAE588", "Vitamins" = "#8DEEEE")
col_cat_nmr = c( "Amino Acid" = "#7AC5CD","Carbohydrate" = "#FFABA9",
                 "Lipids" = "#7CCD7E","Nucleoside" = "#AB82FF", 
                 "Organic Acid" = "#FFA54F", "Other" = "#EEAEEE", 
                 "Peptides" = "#FAE588", "Vitamin" = "#8DEEEE")

#==============================================================================
# Read untargeted data

d.ms <- read.csv("df_combined_POS_NEG.csv")  %>% 
  mutate(Metabolite.name = str_extract(Metabolite.name,"^[^;]+")) 

#Match observed metabolites to HMDB file
d <- read.csv("name_map1.csv")
d2 <- read.csv("name_map2.csv")

d.all <- d %>% 
  filter(!is.na(HMDB)) %>% 
  merge(d2, all = TRUE)

#Load metabolite dataframe
load("HMDB_metabolites.rda")

d.cat <- d.all %>% 
  left_join(metabolites_df %>%
              select(Accession,SuperClass,Class,SubClass),
            by = c("HMDB" = "Accession")) %>% 
  select(c(Query, Match, SuperClass,Class,SubClass ))

#-------------------------------------------------------------------------------
#Average repeat metabolites

#Brute force match repeats
d.ms$Metabolite.name[d.ms$Metabolite.name == "Inosine-5-monophosphate"] <- 
  "Inosine-5'-monophosphate"
d.ms$Metabolite.name[d.ms$Metabolite.name == "2'-Deoxyadenosine-5'-monophosphate"] <-
  "2'-Deoxyadenosine 5'-monophosphate"
d.ms$Metabolite.name[d.ms$Metabolite.name == "(-)-Riboflavin"] <- 
  "Riboflavin"
d.ms$Metabolite.name[d.ms$Metabolite.name == "D-(+)-Pantothenic acid"] <- 
  "Pantothenic acid"

d.ms2 <- d.ms %>% 
  select(c(-Average.Rt.min.,-X, -Adduct.type)) %>% 
  mutate(Metabolite.name = str_to_sentence(Metabolite.name),
         Metabolite = str_extract(Metabolite.name,"^[^;]+") %>% str_trim(),
         Metabolite = ifelse(is.na(str_extract(Metabolite, 
                                               "^[^(]*(?=\\(not validated)")),
                             Metabolite,
                             str_extract(Metabolite,
                                         "^[^(]*(?=\\(not validated)")) 
         %>% str_trim()) %>% 
  pivot_longer(cols = 2:29, names_to = "Hyd", values_to = "val") %>% 
  group_by(Metabolite, Hyd) %>% 
  summarize(avg = mean(val))

#===============================================================================
#Group/organize the categories and subcategories
#Short form/desired groups
#SubClass
Carb = "Carbohydrates and carbohydrate conjugates"
org = "Organic acids and derivatives"

#SuperClass
Lipids = "Lipids and lipid-like molecules"
Nucleosides = "Nucleosides, nucleotides, and analogues"
Other = c("Organic oxygen compounds", "Alkaloids and derivatives","Benzenoids",
          "Organoheterocyclic compounds", "Phenylpropanoids and polyketides")

#Metabolite name
Vitamin = c("4-Pyridoxic acid", "Pyridoxamine", "PYRIDOXINE", "Thiamine",
            "(-)-Riboflavin", "Riboflavin", "NIACIN", "NIACINAMIDE",
            "D-(+)-Pantothenic acid", "Pantothenic acid", "Choline",
            "4-Aminobenzoic acid")
Peptide = c(d.cat$Query[grep("^[A-Za-z]{3}(-[A-Za-z]{3})*$",
                             d.cat$Query)], "Cyclo(Leu-Pro)",
            "Cyclo-prolylglycine", "Isoleucylisoleucine")
AA = c("Tyramine", "Tryptophan", "TRYPTOPHAN", "Arginine", 
       "Asparagine", "Aspartic acid", "Citrulline", "Cystathionine",
       "D-Alloisoleucine", "Glutamate", "Glutamine", "Histidine", 
       "Homoarginine", "Isoleucine", "Isoleucine (not validated)",
       "L-Glutamic acid", "L-Tyrosine", "Ornithine", "Phenylalanine",
       "Serine", "Threonine", "Tyrosine", "Valine", "Carnitine", 
       "5-Hydroxytryptophan", "Kynurenine", "Agmatine", "Histamine",
       "L-Histidinol", "Tryptamine", "Lysine", "Methionine", "Betaine")
Lip = c("FA 18:1+1O", "FA 18:1+2O", "FA 18:1+3O", "Phosphocholine", 
        "1-Decanoyl-2-hydroxy-sn-glycero-3-phosphocholine", 
        "1-Lauroyl-2-hydroxy-sn-glycero-3-phosphocholine")
Nuc = c("Uracil", "Thymine", "Guanine", "Cytosine", "ADENINE",
        "Xanthosine (Not validated)", "Hypoxanthine",
        "Flavin??adenine??dinucleotide")
oth = c("Spermidine", "Spermine", "Secoisolariciresinol")
carbo = c("?-Gentiobiose")

#-------------------------------------------------------------------------------
#Adjust the categories based on manual fixes
d.cat <- d.cat %>% 
  mutate(Category = case_when(
    Query %in% Vitamin ~ "Vitamins",
    Query %in% Peptide ~ "Peptides",
    Query %in% AA ~ "Amino Acid",
    Query %in% Lip ~ "Lipids",
    Query %in% Nuc ~ "Nucleosides",
    Query %in% oth ~ "Other",
    Query %in% carbo ~ "Carbohydrates",
    SubClass %in% Carb ~ "Carbohydrates",
    SuperClass %in% org ~ "Organic acids",
    SuperClass %in% Lipids ~ "Lipids",
    SuperClass %in% Nucleosides ~ "Nucleosides",
    is.na(SuperClass) ~ "Other",
    SuperClass %in% Other ~ "Other", #Take this line out to visualize with more variety
    TRUE ~ SuperClass
  ))

#-------------------------------------------------------------------------------
d.cat.new <- d.cat %>% 
  mutate(Query = str_to_sentence(Query),
         Query = ifelse(is.na(str_extract(Query, 
                                          "^[^(]*(?=\\(not validated)")),
                        Query,
                        str_extract(Query,
                                    "^[^(]*(?=\\(not validated)")))

d.ms.cat <- merge(d.ms2, d.cat.new %>% select(c(Query, Category)),
                  by.x = "Metabolite", by.y = "Query") %>% 
  distinct()

#Visualize overall category breakdown
d.cat.breakdown <- d.ms.cat %>% 
  select(c(Metabolite, Category)) %>% 
  distinct() %>% 
  mutate(Category = as.factor(Category))


ggplot(d.cat.breakdown) +
  geom_bar(aes(x = reorder(Category, Category, function(x) length(x)), fill = Category)) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = -60, hjust = 0, vjust = 1), 
        text = element_text(size = 18), legend.position = "none") +
  scale_fill_manual("Category", values = col_cat)+
  labs(x = "Metabolite Category", y = "Number of Metabolites")

ggsave("Fig1.pdf", width = 6, height = 5)

#===============================================================================
#Compare the categories for the top metabolites within a product
#------------------------------------------------------------------------------
#Product type
hyd.prod <- c(
  HyPep.7504_1 = "Cotton",
  HyPep.7504_2 = "Cotton",
  HyPep.7504_3 = "Cotton",
  HyPea.7404_1 = "Pea",
  HyPea.7404_2 = "Pea",
  HyPea.7404_3 = "Pea",
  HyPep.1510_1 = "Soy",
  HyPep.1510_2 = "Soy",
  HyPep.1510_3 = "Soy",
  HyPep.1510_4 = "Soy",
  HyPep.4601N_1 = "Wheat",
  HyPep.4601N_2 = "Wheat",
  HyPep.4601N_3 = "Wheat",
  HyPep.YE_1 = "YE",
  HyPep.YE_2 = "YE",
  HyPep.YE_3 = "YE",
  HY.YEST.412_1 = "Hy-Yest-412",
  HY.YEST.412_2 = "Hy-Yest-412",
  HY.YEST.412_3 = "Hy-Yest-412",
  HY.YEST.466_1 = "Hy-Yest-466",
  HY.YEST.466_2 = "Hy-Yest-466",
  HY.YEST.466_3 = "Hy-Yest-466",
  HY.YEST.503_1 = "Hy-Yest-503",
  HY.YEST.503_2 = "Hy-Yest-503",
  HY.YEST.503_3 = "Hy-Yest-503",
  HY.YEST.555_1 = "Hy-Yest-555",
  HY.YEST.555_2 = "Hy-Yest-555",
  HY.YEST.555_3 = "Hy-Yest-555"
)

#===============================================================================
#Comparison of the top 50 metabolites
reframe <- d.ms.cat %>% 
  rename(Hydrolysate = Hyd, Val = avg) %>% 
  mutate(Prod = hyd.prod[Hydrolysate],
         Prod = factor(Prod, levels = c("Cotton", "Pea", "Soy", "Wheat", "YE",
                                        "Hy-Yest-412", "Hy-Yest-466",
                                        "Hy-Yest-503", "Hy-Yest-555")))
topmetab <- reframe %>% 
  group_by(Prod, Metabolite) %>% 
  mutate(avg = mean(Val)) %>% 
  ungroup() %>% 
  select(c(Metabolite, Category, Prod, avg)) %>% 
  distinct() %>% 
  group_by(Prod) %>% 
  mutate(Q50 = quantile(avg, probs = .5, na.rm = TRUE),
         Q75 = quantile(avg, probs = .75, na.rm = TRUE)) %>% 
  filter(avg > Q75)
#slice_max(order_by = avg, n = 50)

topmetab <- right_join(reframe, topmetab, 
                       by = c("Metabolite", "Prod", "Category"))

top.prod <- topmetab %>% 
  select(c("Metabolite", "Category", "Prod", "avg")) %>% 
  distinct() %>% 
  mutate(sf = str_sub(Metabolite, 1, 20)) #short form


#Bubble plot of the top metabolites
ggplot(top.prod, aes(x = Prod, y = sf)) +
  geom_count(aes(color = avg, size = avg)) +
  guides(color = guide_legend(title = "Peak Area"), 
         size = guide_legend(title = "Peak Area")) + 
  facet_wrap(~Category, scales = "free")+
  theme_classic() +
  labs(x = "Hydrolysate", y = "Metabolite") +
  theme(#axis.text.y = element_blank(),
    axis.text.x = element_text(angle = -55, hjust = 0, vjust = 1),
    axis.text.y = element_text(angle = 15, hjust = 1, vjust = 0),
    axis.text = element_text(size = 16), 
    axis.title = element_text(size = 20))

#Repeat as separate  plots - selecting the main categories
top.prod.group <- top.prod %>% 
  group_split(Category)

plot.top <- function(d){
  P <- ggplot(d, aes(x = Prod, y = sf)) +
    geom_count(aes(color = avg, size = avg)) +
    guides(color = guide_legend(title = "Peak Area"), 
           size = guide_legend(title = "Peak Area")) + 
    theme_bw() +
    labs(x = "Hydrolysate", y = "Metabolite") +
    theme(axis.text.x = element_text(angle = -55, hjust = 0, vjust = 1),
          axis.text.y = element_text(angle = 15, hjust = 1, vjust = 0),
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14))
}

plots <- lapply(top.prod.group,plot.top)

aa <- plots[[1]]
nuc <- plots[[4]]
pep <- plots[[7]]
vit <- plots[[8]]

plot_grid(aa, nuc, pep, vit, labels = c("A", "B", "C", "D"), label_size = 20)
ggsave("Top_metab.pdf", height = 20, width = 20)

#===============================================================================
#Metabolite heatmap

#-------------------------------------------------------------------------------
#Product labels
prod.expanded <- c(
  HyPep.7504_1 = "Cotton_1",
  HyPep.7504_2 = "Cotton_2",
  HyPep.7504_3 = "Cotton_3",
  HyPea.7404_1 = "Pea_1",
  HyPea.7404_2 = "Pea_2",
  HyPea.7404_3 = "Pea_3",  
  HyPep.1510_1 = "Soy_1",
  HyPep.1510_2 = "Soy_2",
  HyPep.1510_3 = "Soy_3",
  HyPep.1510_4 = "Soy_4",
  HyPep.4601N_1 = "Wheat_1",
  HyPep.4601N_2 = "Wheat_2",
  HyPep.4601N_3 = "Wheat_3",
  HyPep.YE_1 = "YE_1",
  HyPep.YE_2 = "YE_2",
  HyPep.YE_3 = "YE_3",
  HY.YEST.412_1 = "Hy-Yest 412_1",
  HY.YEST.412_2 = "Hy-Yest 412_2",
  HY.YEST.412_3 = "Hy-Yest 412_3",
  HY.YEST.466_1 = "Hy-Yest 466_1",
  HY.YEST.466_2 = "Hy-Yest 466_2",
  HY.YEST.466_3 = "Hy-Yest 466_3",
  HY.YEST.503_1 = "Hy-Yest 503_1",
  HY.YEST.503_2 = "Hy-Yest 503_2",
  HY.YEST.503_3 = "Hy-Yest 503_3",
  HY.YEST.555_1 = "Hy-Yest 555_1",
  HY.YEST.555_2 = "Hy-Yest 555_2",
  HY.YEST.555_3 = "Hy-Yest 555_3"
)

#-------------------------------------------------------------------------------

#Average of the repeat metabolites
d.heat <- d.ms.cat %>% 
  mutate(prod = prod.expanded[Hyd],
         avg = as.numeric(avg),
         prod = factor(prod, 
                       levels = c("Cotton_1", "Cotton_2","Cotton_3","Pea_1",
                                  "Pea_2", "Pea_3", "Soy_1", "Soy_2", "Soy_3", 
                                  "Soy_4", "Wheat_1","Wheat_2", "Wheat_3",
                                  "YE_1", "YE_2", "YE_3", "Hy-Yest 412_1", 
                                  "Hy-Yest 412_2", "Hy-Yest 412_3",
                                  "Hy-Yest 466_1", "Hy-Yest 466_2",
                                  "Hy-Yest 466_3","Hy-Yest 503_1",
                                  "Hy-Yest 503_2", "Hy-Yest 503_3",
                                  "Hy-Yest 555_1","Hy-Yest 555_2",
                                  "Hy-Yest 555_3"))) %>% 
  select(-Hyd) %>%
  arrange(prod) %>% 
  pivot_wider(names_from = "prod", values_from = "avg") %>% 
  as.data.frame()

#Create matrix for mapping
m.int <- as.matrix(d.heat[, -1:-2]) #Remove metabolite column
rownames(m.int) <- c(d.heat[, 1]) #Set row names

# Create row annotations
row_annotation <- data.frame(Category = as.factor(d.heat[,2]))
rownames(row_annotation) <- rownames(m.int)

col_cat2 = list(Category = c( "Amino Acid" = "#7AC5CD","Carbohydrates" = "#FFABA9",
                              "Lipids" = "#7CCD7E","Nucleosides" = "#AB82FF", 
                              "Organic acids" = "#FFA54F", "Other" = "#EEAEEE", 
                              "Peptides" = "#FAE588", "Vitamins" = "#8DEEEE"))

p.heat <- pheatmap(m.int,
                   scale = "row",
                   show_rownames = FALSE, #Improve readability by removing metabolite names
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   color = colorRampPalette(c( "navy","white", "#E31A1C"))(50),
                   annotation_row = row_annotation,
                   angle_col = 315,
                   annotation_colors = col_cat2
)

pdf("metaboliteheat.pdf", height = 15, width = 8)
p.heat
dev.off()

png(file = "metab_heat.png", width = 1200, height = 1400, res = 150)  
p.heat
dev.off()

#===============================================================================
# Batch to batch variability

d.cov <- d.ms2  %>% 
  mutate(Prod = hyd.prod[Hyd]) %>% 
  group_by(Metabolite) %>% 
  mutate(conc_scale = avg/mean(avg)) %>%
  ungroup() %>% 
  group_by(Prod, Metabolite) %>% 
  mutate(batch.cov = sqrt(var(conc_scale))) 

d.cov.quant <- d.cov %>% 
  select(c(Metabolite, Prod, batch.cov)) %>%
  distinct() %>% 
  group_by(Prod) %>% 
  summarize(Q10 = quantile(batch.cov, probs = 0.1, na.rm = TRUE),
            Q25 = quantile(batch.cov, probs = .25, na.rm = TRUE),
            Q50 = quantile(batch.cov, probs = .5, na.rm = TRUE),
            Q75 = quantile(batch.cov, probs = .75, na.rm = TRUE),
            Q90 = quantile(batch.cov, probs = .9, na.rm = TRUE))


d.cov.overall <- d.cov %>% 
  group_by(Prod, Metabolite) %>% 
  mutate(batch.cov = sqrt(var(conc_scale))) %>% 
  ungroup() %>% 
  group_by(Metabolite) %>% 
  mutate(overall = sqrt(var(conc_scale))) %>% 
  select(c(Metabolite,Prod,batch.cov,overall)) %>% 
  distinct() %>% 
  pivot_wider(names_from = Prod, values_from = batch.cov)


#batch metabolite variability in the top 25%
top.sample <- reframe %>% 
  group_by(Hydrolysate) %>% 
  mutate(Q50 = quantile(Val, probs = .5, na.rm = TRUE),
         Q75 = quantile(Val, probs = .75, na.rm = TRUE)) %>% 
  filter(Val > Q75) %>% 
  group_split()

test <- list(top.sample[[1]]$Metabolite, top.sample[[2]]$Metabolite,
             top.sample[[3]]$Metabolite)
venn.diagram(
  x = test,
  category.names = c("A", "B", "C"),
  filename = "Test_yeast.png",
  cex = 1.5, cat.cex = 1.8, margin = 0.1,
  cat.dist = c(0.1, 0.1, 0.08)
)

#===============================================================================
# PCA
m <- pivot_wider(d.ms.cat, names_from = Hyd, values_from = avg)

compounds <- m$Metabolite
samples <- colnames(m)[-1]
category <- m$Category
m <- m[, -1:-2] %>%
  as.matrix() %>%
  t()

colnames(m) <- compounds
pca <- prcomp(m, scale = TRUE)

# Variance covered by PC1 and PC2
var_exp <- summary(pca)$importance[2, ]

d.pca <- pca$x[,1:2] %>% 
  as.data.frame() %>%
  mutate(
    sample = rownames(pca$x),
    Product = hyd.prod[sample]
  )

p.scores <- ggplot(d.pca, aes(x = PC1, y = PC2, colour = Product)) +
  geom_point(size = 5) + 
  theme_classic(20) +
  theme(text = element_text(size = 14), legend.position = "bottom") +
  scale_colour_manual("Sample type", values = Col_all) +
  xlab(paste("PC 1 (", round(var_exp[1], 2)*100, "%)", sep = "")) +
  ylab(paste("PC 2 (", round(var_exp[2], 2)*100, "%)", sep = ""))

d.pca <- pca$rotation[,1:2] %>% 
  as.data.frame() %>%
  mutate(
    compound = rownames(pca$rotation),
    compound_type = category
  )

p.loadings <- ggplot(
  d.pca, aes(label = compound, colour = compound_type)) +
  geom_segment(
    aes(x = 0, y = 0, xend = PC1, yend = PC2), 
    arrow = arrow(length = unit(0.03, "npc"))) +
  geom_text_repel(aes(x = PC1, y = PC2), size = 6) + 
  theme_classic(20) +
  theme(text = element_text(size = 14), legend.position = "bottom") + 
  scale_colour_manual("Category", values = col_cat) +
  xlab(paste("PC 1 (", round(var_exp[1], 2)*100, "%)", sep = "")) +
  ylab(paste("PC 2 (", round(var_exp[2], 2)*100, "%)", sep = ""))


plot_grid(p.scores, p.loadings, labels = c('A', 'B'), label_size = 24, 
          align = "hv")
ggsave("PCA.pdf", width = 16, height = 8)
ggsave("PCA.png", width = 16, height = 8)


