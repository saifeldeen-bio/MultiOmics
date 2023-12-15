library(data.table)
library(MOFA2)
library(ggplot2)
# import data as matrices
setwd("/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/")
circRNA <- "/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/circRNA.csv"
circRNA <- read.csv(circRNA, header = T, sep = ",", row.names = 1)
miRNA <- "/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/miRNA.tsv"
miRNA <- read.csv(miRNA, header = T, sep = "\t", row.names = 1)
phosphorylation <- "/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/phos.tsv"
phosphorylation <- read.csv(phosphorylation, header = T, sep = ",", row.names = 1)
proteome <- "/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/prot.tsv"
proteome <- read.csv(proteome, header = T, sep = ",", row.names = 1)
transcriptome <- "/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/Trans.tsv"
transcriptome <- read.csv(transcriptome, header = T, sep = ",", row.names = 1)
metaData <- "/media/saifeldeen/E43E57ED3E57B6F0/Research_Projects/Private/MultiOmics/Dr_Sally/DATA/merged/meta.csv"  
metaData <- read.csv(metaData, header = T, sep = ",")

circRNA <- as.matrix(circRNA)
miRNA <- as.matrix(miRNA)
phosphorylation <- as.matrix(phosphorylation)
proteome <- as.matrix(proteome)
transcriptome <- as.matrix(transcriptome)


#data <- make_example_data(
#  n_views = 2, 
#  n_samples = 200, 
#  n_features = 1000, 
#  n_factors = 10
#)[[1]]


#sample data
data <- list(circRNA,miRNA, phosphorylation, proteome,transcriptome)
#groups <- metaData$Category
multiOmocs <- create_mofa(data = data)
views_names(multiOmocs) <- c("circRNA","miRNA", "Phosphorylation", "Proteome","Transcriptome")
#if you haave groups (metaData)

# long data frame
#dt = fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/getting_started/data.txt.gz")
#head(dt)
# Let's ignore the grouping information to start
#dt[,group:=NULL]

#Create the MOFA object

#MOFAobject <- create_mofa(dt)
############################### Visualization #####################################
plot_data_overview(multiOmocs)

############################### Define options #######################################
# Define data options
#if groups/views have different ranges/variances, it is good practice to scale each group to unit variance. Default is FALSE
data_opts <- get_default_data_options(multiOmocs)
head(data_opts)
data_opts$scale_views <- TRUE
model_opts <- get_default_model_options(multiOmocs)
model_opts$num_factors <- 10

train_opts <- get_default_training_options(multiOmocs)
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 50000
################################### Build and train the MOFA object ##########################
#Prepare the MOFA object
MOFAobject <- prepare_mofa(
  object = multiOmocs,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
#Train the MOFA model
outfile = file.path(getwd(),"model.hdf5")
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = T)
#environment location: /home/saifeldeen/.cache/R/basilisk/1.10.2/0
# PREFIX=/home/saifeldeen/.cache/R/basilisk/1.10.2/0
################################ Down stream analysis ########################################
#load model
modelPath <- "model.hdf5"
model <- load_model(modelPath)
# overview
plot_data_overview(model)
# Add metadata to the model
#Nsamples = sum(model@dimensions$N)

#sample_metadata <- data.frame(
#  sample = samples_names(model)[[1]],
#  condition = sample(c("A","B"), size = Nsamples, replace = T),
#  age = sample(1:100, size = Nsamples, replace = T)
#)

samples_metadata(model) <- metaData
head(model@samples_metadata, n=3)
slotNames(model) # Slots Name
names(MOFAobject@data)
dim(model@data)
names(model@expectations)
#################################### Plot Factor Correlation #################################
plot_factor_cor(model)
#################################### Variance decomposition ###################################
model@cache[["variance_explained"]] # variance in each view and group and factors
#plot varinace
plot_variance_explained(model, x="view", y="factor")
plot_variance_explained(model, x="group", y="factor", plot_total = T)

#Visualisation of single factors
plot_factor(model, 
            factor = 1:3,
            color_by = "group1"
)

p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "group1",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)
# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red")) +
  scale_fill_manual(values=c("Normal"="blue", "Tumor"="red"))

print(p)

# Visualisation of combinations of factors
plot_factors(model, 
             factors = 1:3,
             color_by = "group1"
)
# Visualisation of feature weights
plot_weights(model,
             view = "Transcriptome",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_top_weights(model,
                 view = "circRNA",
                 factor = 1,
                 nfeatures = 10
)
##################################### Visualisation of patterns in the input data ##########################
# heatmap
plot_data_heatmap(model,
                  view = "Proteome",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)
# Scatter plots
plot_data_scatter(model,
                  view = "Proteome",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "group1"
)
# Non-linear dimensionality reduction
#Run UMAP and t-SNE
set.seed(42)

umap <- run_umap(model)
tsne <- run_tsne(model)
plot_dimred(model,
            method = "TSNE",  # method can be either "TSNE" or "UMAP"
            color_by = "condition"
)
################################ Other functionalities ##########################
#Renaming dimensions
views_names(model) <- c("Transcriptomics", "Proteomics")
factors_names(model) <- paste("Factor", 1:model@dimensions[["K"]], sep=" ")
#Extracting data for downstream analysis
# "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
#Extract factors
factors <- get_factors(model, factors = "all")
lapply(factors,dim)

#Extract weights
# "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(model, views = "all", factors = "all")
lapply(weights,dim)
############################ Enrichment Analysis #####################################
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFA2)
library(MOFAdata)
#ReactomeGS (human):
data("reactomeGS")
head(rownames(reactomeGS), n=3)
#MSigDB 6.0 (human):
# C2: curated gene sets from online pathway databases, publications in PubMed, and knowledge of domain experts.
data("MSigDB_v6.0_C2_human") 

# C5: extracted from the Gene Ontology data.base
data("MSigDB_v6.0_C5_human") 

head(rownames(MSigDB_v6.0_C2_human), n=3)
#MSigDB 6.0 (mouse):

# C2: curated gene sets from online pathway databases, publications in PubMed, and knowledge of domain experts.
data("MSigDB_v6.0_C2_mouse") 

# C5: extracted from the Gene Ontology data.base
data("MSigDB_v6.0_C5_mouse") 

head(rownames(MSigDB_v6.0_C2_mouse), n=3)

model <- load_model("model.hdf5")
features_names(model)[["view_0"]] <- toupper(features_names(model)[["view_0"]])
head(features_names(model)[["view_0"]])

#runGSEA
enrichment.parametric <- run_enrichment(model,
                                        view = "view_0", factors = 1:3,
                                        feature.sets = MSigDB_v6.0_C5_mouse,
                                        sign = "negative",
                                        statistical.test = "parametric"
)