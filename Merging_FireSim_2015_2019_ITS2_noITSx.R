### This is the R code associated with merging the 1-year and 5-year post fire
### phyloseq objects for the paper
### Resilience not yet apparent in soil fungal communities of the 
### boreal forest from one to five years after wildfire across a severity gradient###

# Paper authors: Thea Whitman, Jamie Woolet, Miranda C. Sikora, Dana B. Johnson, Denyse A. Dawe, and Ellen Whitman 

# © 2025. This work is openly licensed via CC BY NC SA 4.0.
# https://creativecommons.org/licenses/by-nc-sa/4.0/

# Load packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(DECIPHER)
library(Biostrings)

# Import 2019 dataset
ps.5 = readRDS("WB2019.ITS2.noITSx.cutadapt.ps")
ps.5

# Bringing in the 1 year post fire dataset
ps.1 = readRDS("WB2015.ITS2.noITSx.cutadapt.ps")
ps.1

# All 2019 samples are in the 2015 dataset; some 2015 samples not in 2019
sample_data(ps.1)$Site_ID %in% sample_data(ps.5)$Site_ID
sample_data(ps.5)$Site_ID %in% sample_data(ps.1)$Site_ID

# Keep only sites that are in the 2019 dataset
ps.1 = subset_samples(ps.1,Site_ID %in% sample_data(ps.5)$Site_ID)

# Check for same samples in both datasets
sample_data(ps.1)$Sample_ID %in% sample_data(ps.5)$Sample_ID
sample_data(ps.5)$Sample_ID %in% sample_data(ps.1)$Sample_ID
sample_data(ps.5)$Sample_ID[30]
# Sample 15S-NT-U06M likely was not captured in 2015, possibly just due to spatial variability
# It's a treed wetland, so O horizon may have been too thick.
# We'll drop it from the paired dataset.
ps.5 =  subset_samples(ps.5, !(Sample_ID %in% c("15S-NT-U06M")))

# Note: WB19-38-B was the new site established because the previous site had been logged between sampling years. (15S-NT-38)
# If we check the ordinations, the M and O horizons for the B and non-B sites are similar to each other.
# We'll use the B site going forward, since the first site has likely been disturbed
# We also have WB19-27-O-1 and WB19-27-O-2.
# It seems like WB19-27-O-2 is the correct sample for that site.
ps.5 =  subset_samples(ps.5, !(sample_names(ps.5) %in% c("15S-NT-27-O-1","15S-NT-38-M","15S-NT-38-O")))

# Remove two unburned mineral samples that we only collected in 2015 (but do have O horizons for)
ps.1 =  subset_samples(ps.1, !(Sample_ID %in% c("15S-NT-U07M","15S-NT-U10M")))

# Make sure no zero-count samples are present
ps.5 = subset_taxa(ps.5, taxa_sums(ps.5)>0)
ps.1 = subset_taxa(ps.1, taxa_sums(ps.1)>0)

# Check phyloseq objects
ps.5
ps.1

# We do seem to have more OTUs in the 2019 dataset. That could indicate recovery, but it could also be due to seq depth or small
# changes in sampling depth for mineral horizon between years. Since O horizon sampling stayed the same, if that were the case,
# we should not see diffs in O horizon.
ps.5.O = subset_taxa(ps.5, sample_data(ps.5)$Org_or_Min=="O")
ps.1.O = subset_taxa(ps.1, sample_data(ps.1)$Org_or_Min=="O")
ps.5.O = subset_taxa(ps.5.O, taxa_sums(ps.5.O)>0)
ps.1.O = subset_taxa(ps.1.O, taxa_sums(ps.1.O)>0)
ps.5.O
ps.1.O

# No, it persists in O horizon - maybe even more so. 
# Makes sense that fungi might be highly affected in O horizon, 1 year post fire.
# May suggest recovery.

# Add year burned
sample_data(ps.1)$Years_Since_Fire = "1"
sample_data(ps.5)$Years_Since_Fire = "5"

# Probably do want to throw out any taxa that didn't get a phylum match
# Might be overly conservative, but since we didn't use ITSx, it's best
# to avoid accidentally including other eukaryotes.

sum(is.na(tax_table(ps.5)[,"Phylum"])) / ntaxa(ps.5)
sum(is.na(tax_table(ps.1)[,"Phylum"])) / ntaxa(ps.1)
ps.5 = subset_taxa(ps.5, !is.na(tax_table(ps.5)[,"Phylum"]))
ps.1 = subset_taxa(ps.1, !is.na(tax_table(ps.1)[,"Phylum"]))


############## First approach is to take the fasta files and merge based on identical sequences

# Because the sequenced regions were identical between the two years,
# we should be able merge all OTUs by their sequence alone

# Bring in the DNA sequences
dna.5 = readDNAStringSet("../Seq-processing/sequences.2019.noITSx.cutadapt.fasta")
dna.1 = readDNAStringSet("../Seq-processing/sequences.2015.noITSx.cutadapt.fasta")

# Extract DNA sequences
seqString.5 = c()
for (i in 1:length(dna.5)){
  seq = paste(dna.5[[i]])
  seqString.5[i]=seq
}
head(seqString.5)
length(seqString.5)

# Make table with ID and fasta
df.5 = data.frame(OTU=dna.5@ranges@NAMES,Sequence=seqString.5)
dim(df.5)
head(df.5)

# Subset it so only OTUs that remain in the phyloseq object are present
df.5 = df.5[df.5$OTU %in% taxa_names(ps.5),]
dim(df.5)

# Let's rename our OTUs.
ps.5.fasta = ps.5
sum(taxa_names(ps.5.fasta)==df.5$OTU)

# They are in the same order, so we can do the following
taxa_names(ps.5.fasta)=df.5$Sequence

### Follow the same procedure for 2015 dataset

# Extract DNA sequences
seqString.1 = c()
for (i in 1:length(dna.1)){
  seq = paste(dna.1[[i]])
  seqString.1[i]=seq
}
length(seqString.1)
# Have a string of the DNA sequences

# Make table with ID and fasta
df.1 = data.frame(OTU=dna.1@ranges@NAMES,Sequence=seqString.1)
dim(df.1)
head(df.1)

# Subset it so only OTUs that remain in the phyloseq object are present
df.1 = df.1[df.1$OTU %in% taxa_names(ps.1),]
dim(df.1)

# Let's rename our OTUs.
ps.1.fasta = ps.1
sum(taxa_names(ps.1.fasta)==df.1$OTU)

# They are in the same order, so we can do the following
taxa_names(ps.1.fasta)=df.1$Sequence

### Now we have two phyloseq objects where the names of the OTUs are just their DNA sequences.

# Ok, we now should have taxa names are the fasta sequence for each phyloseq object.
# Now we want to merge them.
ps.merged.fasta = merge_phyloseq(ps.1.fasta,ps.5.fasta)
ps.merged.fasta
ps.1.fasta
ps.5.fasta
# All the samples made it through.

# How many OTUs were there in the two datasets alone?
ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)

# After merging, how many OTUs were there?
ntaxa(ps.merged.fasta)

# What was the difference?
ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)-ntaxa(ps.merged.fasta)

# This is a relatively small fraction of the total OTUs - 19%
((ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)-ntaxa(ps.merged.fasta))) /ntaxa(ps.merged.fasta)

# If we compare it just to the 2015 dataset, it's higher (56% of OTUs)
# This is much better than the initial bioinformatics with ITSx
((ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)-ntaxa(ps.merged.fasta))) /ntaxa(ps.1)

# We can check fraction of total sequences it represents
ps.1.fasta.norm = transform_sample_counts(ps.1.fasta, function(x) x/sum(x))
mdf.ps.1.fasta.norm = psmelt(ps.1.fasta.norm)
shared.taxa = taxa_names(ps.1.fasta)[taxa_names(ps.1.fasta) %in% taxa_names(ps.5.fasta)]

df = mdf.ps.1.fasta.norm%>%
  filter(OTU %in% shared.taxa)%>%
  group_by(Sample)%>%
  summarize(Sum=sum(Abundance))
df

hist(df$Sum)
mean(df$Sum)
# Previously, even accounting for abundance, there was low representation (13%)
# Now, we have 89% of reads shared

ps.5.fasta.norm = transform_sample_counts(ps.5.fasta, function(x) x/sum(x))
mdf.ps.5.fasta.norm = psmelt(ps.5.fasta.norm)
shared.taxa = taxa_names(ps.5.fasta)[taxa_names(ps.5.fasta) %in% taxa_names(ps.1.fasta)]

df = mdf.ps.5.fasta.norm%>%
  filter(OTU %in% shared.taxa)%>%
  group_by(Sample)%>%
  summarize(Sum=sum(Abundance))
df

hist(df$Sum)
mean(df$Sum)
# For the 2019 dataset, it's 67% of all reads - certainly less than the 2015 dataset.

hist(sample_sums(ps.5.fasta))
hist(sample_sums(ps.1.fasta))
mean(sample_sums(ps.5.fasta))
mean(sample_sums(ps.1.fasta))
sum(sample_sums(ps.5.fasta))
sum(sample_sums(ps.1.fasta))
sum(sample_sums(ps.1.fasta))/sum(sample_sums(ps.5.fasta))
# The 2019 library was sequenced more deeply than the 2015 library,
# even though each sample was extracted and sequenced twice.
# This is another potential contributor to the increase in observed OTUs in 2019.

# At a minimum, we might check to merge sequences that are 100%ID, but just different lengths,
# in case there were small variations in trimming

################## Another approach to merging ###########################


# Because the ITS region is of variable lengths,
# and in case the trimming from ITSx may have trimmed at slightly different locations,
# we merge all OTUs that are identical except for length.

# Running clustering on full fasta file
# Set processors

# Collect the DNA sequences
# We can use the OTUs in the merged object created above
dna = DNAStringSet(taxa_names(ps.merged.fasta))

# Align the sequences
# Since we're clustering at 0%, the exact alignment isn't super critical.
aln = DECIPHER::AlignSeqs(dna, processors=NULL, iterations=2, refinements=1)
#saveRDS(aln,"ps.merged.fasta.alignment")

#DECIPHER::BrowseSeqs(aln, highlight=0)

# Find similarity between all seqs
d = DECIPHER::DistanceMatrix(aln, processors = NULL, includeTerminalGaps = FALSE,
                             penalizeGapLetterMatches = TRUE,
                             penalizeGapGapMatches = FALSE)

# Cluster at identical similarity
clusters = DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.0, # must be identical (except ends); could use `cutoff = 0.03` for a 97% OTU 
  processors = NULL)
head(clusters)
length(unique(as.factor(clusters$cluster)))
dim(clusters)

# This would bring us down to 4138 OTUs instead of 4841 OTUs
# That's a decent reduction.
# Let's see how common the clusters are
clust.num = clusters%>%
  dplyr::group_by(cluster)%>%
  dplyr::summarize(TotalASVs=n())%>%
  dplyr::arrange(-TotalASVs)
clust.num[1:10,]

# The max number of OTUs merged at 100% ID was 1
# So, no OTUs had identical sequences but different lengths.
# With previous appraoches, there were a number of these - this is interesting with so many OTUs.
# Suggests this approach might be better than previous ones - less variable trimming

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters = clusters %>%
  dplyr::mutate(OTUID = row.names(clusters))

head(clusters)
# The DNA seqs were named by numbers
# We need to get the names back from the phyloseq object
# In the ps object, they are just the sequences. We should be able to extract from the dna object.
clusters$taxanames = names(dna)

# The OTU names for the cluster IDs are in the same order as the OTU names in the phyloseq object. Good.
sum(taxa_names(ps.merged.fasta) == clusters$taxanames)

# Get the taxonomy table and add the clustering data
tax.tab.clust = as.matrix(tax_table(ps.merged.fasta))
tax.tab.clust = cbind(tax.tab.clust,as.matrix(clusters))
head(tax.tab.clust)
tax.tab.clust = tax_table(tax.tab.clust)
tax_table(ps.merged.fasta) = tax.tab.clust
tax_table(ps.merged.fasta)[,"cluster"]

ps.glommed = tax_glom(ps.merged.fasta,"cluster")
ps.glommed
length(unique(tax_table(ps.glommed)[,"cluster"]))
# It looks like we have successfully merged all our clusters
# Not really needed as there was no clustering in the end

# We could additionally merge by species
ps.glommed.spp = tax_glom(ps.glommed,"Species",NArm=FALSE)
ps.glommed.spp
ps.glommed

# This loses the cluster ID.
tax = data.frame(tax_table(ps.glommed.spp))
tax$OTU = row.names(tax)
Clusters = data.frame(tax_table(ps.glommed)[,"cluster"])
Clusters$OTU = row.names(Clusters)
tax = merge(tax,Clusters,by="OTU",all.x=TRUE)
tax = tax[,c("Kingdom","Phylum","Class","Order","Family","Genus","Species","cluster.y","OTU")]
colnames(tax)=c("Kingdom","Phylum","Class","Order","Family","Genus","Species","cluster","OTU")
tax2 = tax_table(tax)
taxa_names(tax2)=tax$OTU
tax2 = tax2[,1:8]
colnames(tax2)=c("Kingdom","Phylum","Class","Order","Family","Genus","Species","cluster")
head(tax2)
# Add back into glommed
tax_table(ps.glommed.spp)=tax2

# Including NArm=FALSE keeps the taxa with no species assignment in, but it still groups them
# So, we want to remove all these taxa from the dataset, and replace them with the original OTUs.
NASpeciesOTUs = taxa_names(subset_taxa(ps.glommed.spp,is.na(tax_table(ps.glommed.spp)[,"Species"])))
ps.glommed.spp.noNAs = subset_taxa(ps.glommed.spp,!(taxa_names(ps.glommed.spp) %in% NASpeciesOTUs))
ps.glommed.spp.noNAs

# This gives us the glommed dataset without the NA species taxa
# Now we want to get all the taxa from the previous phyloseq object that did have NA in species.
NASpeciesOTUs2 = taxa_names(subset_taxa(ps.glommed,is.na(tax_table(ps.glommed)[,"Species"])))

ps.merged.allNAs = subset_taxa(ps.glommed,(taxa_names(ps.glommed) %in% NASpeciesOTUs2))
ps.merged.allNAs
tax_table(ps.merged.allNAs) = tax_table(ps.merged.allNAs)[,1:8]
length(NASpeciesOTUs2)

# Merge the two
newtaxtab = tax_table(rbind(data.frame(tax_table(ps.glommed.spp.noNAs)),data.frame(tax_table(ps.merged.allNAs))))
rownames = c(taxa_names(ps.glommed.spp.noNAs),taxa_names(ps.merged.allNAs))
rownames(newtaxtab) = rownames
colnames(newtaxtab) = colnames(tax_table(ps.glommed.spp.noNAs))

newotutab = otu_table(rbind(data.frame(t(otu_table(ps.glommed.spp.noNAs))),data.frame(t(otu_table(ps.merged.allNAs)))),taxa_are_rows=TRUE)

newsamdat = sample_data(ps.glommed.spp.noNAs)

head(sample_names(newsamdat))
head(sample_names(newotutab))
tail(sample_names(newsamdat))
tail(sample_names(newotutab))

sample_names(newotutab)=sample_names(newsamdat)

ps.glommed.spp.wNAs = phyloseq(newtaxtab,newotutab,newsamdat)

ps.glommed.spp
ps.glommed
ps.glommed.spp.wNAs
# We added back the OTUs without species IDs that we'd lost

# Now we have an OTU table where taxa that are a 100% ID match are joined as the same OTU
# and taxa that were given the same species name are also joined as the same OTU.
# There is no clustering at this point, other than ignoring insertions or deletions at the 
# end of the sequence

# Let's see the fraction of shared reads now.

ps.5.aligned = subset_samples(ps.glommed.spp.wNAs,sample_data(ps.glommed.spp.wNAs)$Years_Since_Fire=="5")
ps.5.aligned = subset_taxa(ps.5.aligned, taxa_sums(ps.5.aligned)>0)
ps.5.aligned
ps.1.aligned = subset_samples(ps.glommed.spp.wNAs,sample_data(ps.glommed.spp.wNAs)$Years_Since_Fire=="1")
ps.1.aligned = subset_taxa(ps.1.aligned, taxa_sums(ps.1.aligned)>0)
ps.1.aligned

ps.5.aligned.norm = transform_sample_counts(ps.5.aligned, function(x) x/sum(x))
mdf.ps.5.aligned.norm = psmelt(ps.5.aligned.norm)
shared.taxa = taxa_names(ps.5.aligned)[taxa_names(ps.5.aligned) %in% taxa_names(ps.1.aligned)]

df = mdf.ps.5.aligned.norm%>%
  filter(OTU %in% shared.taxa)%>%
  group_by(Sample)%>%
  summarize(Sum=sum(Abundance))
df

hist(df$Sum)
mean(df$Sum)

ps.1.aligned.norm = transform_sample_counts(ps.1.aligned, function(x) x/sum(x))
mdf.ps.1.aligned.norm = psmelt(ps.1.aligned.norm)
shared.taxa = taxa_names(ps.1.aligned)[taxa_names(ps.1.aligned) %in% taxa_names(ps.5.aligned)]

df = mdf.ps.1.aligned.norm%>%
  filter(OTU %in% shared.taxa)%>%
  group_by(Sample)%>%
  summarize(Sum=sum(Abundance))
df

hist(df$Sum)
mean(df$Sum)

# Now we have 96% of reads observed in 2015 are from OTUs also present in 2019; 
# 80% of reads from 2019 were also present in 2015.

# Seems basically reasonable.

#saveRDS(ps.glommed.spp.wNAs,"ps.merged.aligned.0pct.glommed.spp.wNAs.noITSx.cutadapt")
