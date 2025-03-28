### This is the R code associated with creating the 1-year post fire
### phyloseq object for the paper
### Resilience not yet apparent in soil fungal communities of the 
### boreal forest from one to five years after wildfire across a severity gradient###

# Paper authors: Thea Whitman, Jamie Woolet, Miranda C. Sikora, Dana B. Johnson, Denyse A. Dawe, and Ellen Whitman 

# © 2025. This work is openly licensed via CC BY NC SA 4.0.
# https://creativecommons.org/licenses/by-nc-sa/4.0/

library(qiime2R)
library(phyloseq)
library(boxrdrive)
library(dplyr)

# Importing each element
otus = qza_to_phyloseq("../Seq-processing/table.2015.noITSx.cutadapt.qza")
samdat = read.csv("../data/WB2019-2015-merged-metadata.csv", header=TRUE)
taxtab = read_qza("../Seq-processing/taxonomy.2015.noITSx.cutadapt.qza")
taxtab = parse_taxonomy(taxtab$data)
repseqs = read_qza("../Seq-processing/rep-seqs.2015.noITSx.cutadapt.qza")

# Creating otu table
otus = otu_table(otus,taxa_are_rows = TRUE)

# Creating sample data
samdat = sample_data(samdat)
sample_names(samdat)=samdat$X.SampleID

# Creating taxonomy table
taxtab2 = tax_table(taxtab)
taxa_names(taxtab2)=row.names(taxtab)
colnames(taxtab2) = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

# Checking compatibility
sample_names(otus) %in% sample_names(samdat)

# Creating phyloseq object
ps = phyloseq(otus,samdat,taxtab2)
ps

# Remove blanks
ps=prune_samples(!sample_data(ps)$Sample_Name=="Blank",ps)
ps

# Remove low-abundance samples
sort(sample_sums(ps))
cutoff = 1000
ps = subset_samples(ps,sample_sums(ps) >cutoff)
ps

# Remove open wetland samples (not present in 2019 dataset)
ps = subset_samples(ps,sample_data(ps)$Veg_Comm != "Open Wetland")
ps

# Merge lab replicate DNA extractions (-1 and -2)
# Note:
# Samples designated A and B are field replicates from the same site
# Sites 20 and 31 are wetlands, so won't be included
# Site 9 R is the correct sample for site 9
ps.merged = merge_samples(ps,sample_data(ps)$Sample_ID)
ps.merged

hist(sample_sums(ps.merged))
min(sample_sums(ps.merged))
# Looks fine

# Need to add sample metadata back
samdat = sample_data(ps)

# Get only one copy for each sample
samdat2 = data.frame(samdat) %>%
  filter(Lab_Rep =="2")%>%
  select(!"X.SampleID")%>%
  select(!"Sample_Name")
c(samdat2$Sample_ID)

# Check that all samples are still present
samdat2$Sample_ID %in% sample_names(sample_data(ps.merged))
sample_names(sample_data(ps.merged)) %in% samdat2$Sample_ID 
sample_names(sample_data(ps.merged)) == samdat2$Sample_ID 

# However, not in same order.
samdat3 = data.frame(sample_data(ps.merged))
samdat3 = merge(samdat3,samdat2,by="Sample_ID",all.y=TRUE)

head(samdat3)
dim(samdat3)
dim(sample_data(ps.merged))

sample_names(sample_data(ps.merged)) == samdat3$Sample_ID 
colnames(samdat3)
samdat3 = samdat3[,c(4,1,75:143)]
colnames(samdat2)
colnames(samdat3)
colnames(samdat3) = colnames(samdat2)
row.names(samdat3) = samdat3$Sample_ID
sample_data(ps.merged) = sample_data(samdat3)

# Should now have samples, merged by lab reps by summing in OTU table,
# with full sample data.

# Formatting variables
sample_data(ps.merged)$Years_Since_Burn = as.factor(sample_data(ps.merged)$Years_Since_Burn )
sample_data(ps.merged)$Severity_Class = factor(sample_data(ps.merged)$Severity_Class,
                                                                         levels = c("High", "Moderate", "Low", 'Unburned'), ordered = TRUE)

# This should be basically ready to go.
ps.merged = prune_taxa(taxa_sums(ps.merged)>0,ps.merged)
ps.merged
ps
# Dropped a few taxa with the removal of open wetlands, which aren't in the 2019 dataset.

# Checking taxa with assignments
sum(is.na(tax_table(ps.merged)[,"Family"]))
sum(!is.na(tax_table(ps.merged)[,"Family"]))
sum(is.na(tax_table(ps.merged)[,"Genus"]))
sum(!is.na(tax_table(ps.merged)[,"Genus"]))
sum(is.na(tax_table(ps.merged)[,"Species"]))
sum(!is.na(tax_table(ps.merged)[,"Species"]))
sum(is.na(tax_table(ps.merged)[,"Phylum"]))

#saveRDS(ps.merged,"WB2015.ITS2.noITSx.cutadapt.ps")
