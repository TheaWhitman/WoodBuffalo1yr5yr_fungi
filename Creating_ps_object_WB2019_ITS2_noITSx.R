### This is the R code associated with creating the 5-year post fire
### phyloseq object for the paper
### Resilience not yet apparent in soil fungal communities of the 
### boreal forest from one to five years after wildfire across a severity gradient###

# Paper authors: Thea Whitman, Jamie Woolet, Miranda C. Sikora, Dana B. Johnson, Denyse A. Dawe, and Ellen Whitman 

# © 2025. This work is openly licensed via CC BY NC SA 4.0.
# https://creativecommons.org/licenses/by-nc-sa/4.0/

library(qiime2R)
library(phyloseq)
library(boxrdrive)
library(Biostrings)

# Importing each element
otus = qza_to_phyloseq("../Seq-processing/table.2019.noITSx.cutadapt.qza")
samdat = read.csv("../data/WB2019-metadata.csv", header=TRUE)
taxtab = read_qza("../Seq-processing/taxonomy.2019.noITSx.cutadapt.qza")
repseqs = read_qza("../Seq-processing/rep-seqs.2019.noITSx.cutadapt.qza")
taxtab = parse_taxonomy(taxtab$data)

# No longer need to add filtering step to remove OTUs not compatible with 2015 dataset
# which required more filtering due to quality profile
# because similar lengths after trimming
# test = DNAStringSet(repseqs$data)
# cutoff = 410
# ShortOTUs = test[width(test)<cutoff]
# ShortOTUs=names(ShortOTUs)

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

# Check out some of the no-phylum taxa
sum(is.na(tax_table(ps)[,"Phylum"]))
Phylum = tax_table(ps)[,"Phylum"]
head(Phylum[is.na(Phylum)])

DNA = repseqs$data
as.character(DNA["aeb4a00dc81be55f1980ba25a0460f57"])
# Generally seem to have relatively poor blast matches in NCBI, although to fungi.
# Could see how influential removing them is on the results

# Remove blanks
ps=prune_samples(!sample_data(ps)$Sample_Name=="Blank",ps)
ps


# Set severity classes as ordered factor
sample_data(ps)$Severity_Class = sample_data(ps)$Severity_Class = factor(sample_data(ps)$Severity_Class,
                                                                         levels = c("High", "Moderate", "Low", 'Unburned'), ordered = TRUE)
# Sample data are initially identical to the 2015 values
# We re-measured pH and C and N for the 2019 samples
# Need to replace those in metadata file

# Importing new pH values for 2019 soils
pH.2019 = read.csv(box_drive("WhitmanLabMaster/WhitmanLab/Projects/WoodBuffalo/NWT-Fire-Microbes/WB2019/data/Soils_data/Summaries/WB2019 pH means_04-30-21.csv"))
# Make sure same sample IDs are present
as.factor(pH.2019$Sample_ID) %in% sample_data(ps)$Sample_ID
sample_data(ps)$Sample_ID%in% as.factor(pH.2019$Sample_ID)
# Yes, good, all sample IDs in our phyloseq object are found in the pH table
# we have a few measured pH values that aren't in our dataset - that's fine

# Check out the ones in the pH table that are not in the phyloseq object and vice versa
pH.2019[!(as.factor(pH.2019$Sample_ID) %in% sample_data(ps)$Sample_ID),]

# 15S-NT-29O isn't in the 5-year sequencing dataset - site 29
# Note: WB19-38-B was the new site established because the previous site had been logged between sampling years. (15S-NT-38)
# So, for the purposes of this paper, 38 might still be better comparison.

# Updating pH data
add.pH = data.frame(sample_data(ps)[(sample_data(ps)$Sample_ID %in% pH.2019$Sample_ID),])
add.pH = merge(add.pH,pH.2019,by="Sample_ID")
colnames(add.pH)
row.names(add.pH) = add.pH$X.SampleID
add.pH$pH = add.pH$Mean_pH
add.pH = add.pH[,1:70]
samdat = sample_data(add.pH)
sample_data(ps)=samdat
sample_data(ps)$pH

# Similarly, update soil C values
C.2019 = read.csv(box_drive("WhitmanLabMaster/WhitmanLab/Projects/WoodBuffalo/NWT-Fire-Microbes/WB2019/data/Soils_data/Summaries/WB2019_MandOhorizon_TCandTN_Summary_06-16-2021.csv"))

# Make sure same sample IDs are present
as.factor(C.2019$Sample_ID) %in% sample_data(ps)$Sample_ID
sample_data(ps)$Sample_ID %in% as.factor(C.2019$Sample_ID)
# As with pH, all sample IDs in our phyloseq object are found in the C table

# Same samples in the CN table that are not in the phyloseq object
C.2019[!(as.factor(C.2019$Sample_ID) %in% sample_data(ps)$Sample_ID),]

# Updating C data
add.C = data.frame(sample_data(ps)[(sample_data(ps)$Sample_ID %in% C.2019$Sample_ID),])
add.C = merge(add.C,C.2019,by="Sample_ID")
colnames(add.C)
row.names(add.C) = add.C$X.SampleID
add.C$TC_pct = add.C$TC_dry_pct
add.C$Total_N_pct = add.C$TN_dry_pct
# Keep just the original columns
add.C = add.C[,1:70]
samdat = sample_data(add.C)
sample_data(ps)=samdat

# Save phyloseq object
#saveRDS(ps,"WB2019.ITS2.noITSx.cutadapt.ps")

