---
  title: "Microbiome Analysis Demo"
author: "Anthony Bonacolta"
date: "`r Sys.Date()`"
output: html_document
---
  # This R Markdown includes code chunks for processing Hailey and Ella's Highschool Coral Microbiome Project
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_knit$set(root.dir = '/Users/ellawilson/Desktop/trimmed') #Set this to the folder on your computer where you will keep everything (raw data + analysis)
```

# 18S DADA2
- This first step (performed by Anthony) trims the raw reads of sequencing adapters and our primers using a program called cutadapt
```{bash}
cd /Users/anthony-bonacolta/Desktop/ColoradoHS_CoralMicrobiome/data

ls *.raw_1.fastq.gz | cut -f 1-2 -d "." > samples
for sample in $(cat samples)
do
echo "On sample: $sample"
cutadapt -a NNNNNNNNCYGCGGTAATTCCAGCTC...CRAAGAYGATYAGATACCRT \
-A NNNNNNNNAYGGTATCTRATCRTCTTYG...GAGCTGGAATTACCGCRG \
--match-read-wildcards \
-m 150:150 -M 500:500 --discard-untrimmed \
-o ${sample}_L001_R1_001_trimmed.fastq.gz -p ${sample}_L001_R2_001_trimmed.fastq.gz \
${sample}.raw_1.fastq.gz ${sample}.raw_2.fastq.gz \
>> cutadapt_primer_trimming_stats.txt 2>&1
done

paste samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")
```

- Now you can process the rest of the workflow in R
- The below code finds the reads in your computer and loads them into R
```{r}
library(dada2); packageVersion("dada2") # This loads DADA2 which we use for all the below functions.
path <- "/Users/ellawilson/Desktop/trimmed" # CHANGE ME to the directory containing the trimmed fastq files
list.files(path)
fnFs <- sort(list.files(path, pattern="_L001_R1_001_trimmed.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001_trimmed.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)
```

- Now we will check the quality of the reads. From these plots we want to make sure our reads. In gray-scale is a heat map of the frequency of each quality score at each base position. The mean quality score at each position is shown by the green line, and the quartiles of the quality score distribution by the orange lines. The red line shows the scaled proportion of reads that extend to at least that position (this is more useful for other sequencing technologies, as Illumina reads are typically all the same length, hence the flat red line).
- Reverse reads tend to be worse than forward reads. 
- The idea here is to find a good trunclen (where we cut off the reads before merging) for the filterandTrim function which comes n ext. This tends to be where the quality of the reads drops significantly in the plots. 
```{r}
plotQualityProfile(fnFs[1:3])
plotQualityProfile(fnRs[1:3])
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```
- Now using the info above we fill in truncLen=c(FXXX,RXXX), with where we see a drop off in our reads. So if it drops off at 220 on the forwards and 190 in the reverses, we'd put truncLen=c(220,190),
```{r}
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=0, truncLen=c(225,210),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE, matchIDs=TRUE)
head(out)
```
- Check the above output to ensure we did not lose too much data with our trimming. We ideally want to keep over 70% of our data here. 

- Now we proceed with the dada2 workflow with these mostly consistent parameters that don't need changing
```{r}
set.seed(100)
errF <- learnErrors(filtFs, nbases=1e8, multithread=2, randomize=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=2, randomize=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names
ddFs <- dada(derepFs, err=NULL, selfConsist=TRUE)
ddRs <- dada(derepRs, err=NULL, selfConsist=TRUE)
plotErrors(ddFs)
plotErrors(ddRs)
dadaFs <- dada(derepFs, err=ddFs[[1]]$err_out, pool=TRUE, multithread=2)
dadaRs <- dada(derepRs, err=ddRs[[1]]$err_out, pool=TRUE, multithread=2)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])
seqtab.all <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab.all)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths.
```
- We now have a table to ASVs (Amplicon Sequence Variants), consider these as essentially "species" of microbes within your dataset.
- Now we want to annotate our ASVs using a database of protists sequences to match them with what they most likely are. 
- In this case we use PR2 (The protist ribosomal rna database)- https://pr2-database.org/
  - Download the latest dada2 compatible file here: https://github.com/pr2database/pr2database/releases 
- Note that this step takes a while (12-48 hours sometimes, so keep your computer open and charged)
  ```{r}
ref_fasta <- "/Users/ellawilson/Desktop/trimmed/pr2_version_5.1.0_SSU_dada2.fasta.gz" #This is downloaded from online, check the microbiome analysis workflow pdf for the current link
taxa <- assignTaxonomy(seqtab, refFasta=ref_fasta, multithread=2, taxLevels = c("Domain","Supergroup","Division", "Subdivision", "Class", "Order","Family","Genus","Species")) # Note the different taxonomic levels here as compared to the 16S workflow
colnames(taxa) <- c("Domain","Supergroup","Division", "Subdivision", "Class", "Order","Family","Genus","Species")
taxa.print <- taxa
```

- Now we output all of dada2 data into excel files which we can read and import into Phyloseq for further analysis! 
  - Adjust the output path accordingly.
```{r}
asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "/PATH/TO/WORKING_DIRECTORY/ASVs.fa")

asv_tab <- t(seqtab)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "/PATH/TO/WORKING_DIRECTORY/ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "/PATH/TO/WORKING_DIRECTORY/ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

asv_tabtax <- cbind(asv_tab,asv_tax)
row.names(asv_tabtax) <- sub(">", "", asv_headers)
write.table(asv_tabtax, "/PATH/TO/WORKING_DIRECTORY/ASVs_counts_taxonomy.tsv", sep="\t", quote=F, col.names=NA)
```
- Great, the read processing portion is done and now you can work with the output files for data vizualization

# Phyloseq 
- Phyloseq is an amazing program which allows us to process and vizualize our microbiome data output from DADA2
- This first step involves importing our data back into R for phyloseq to read
- Make sure you also have your metadata file, telling phyloseq what each sample is.
```{r}
library(phyloseq)
SV_18S <- read.table("/PATH/TO/WORKING_DIRECTORY/ASVs_counts.tsv", row.names = 1, check.names= FALSE, sep = "\t", header = TRUE)
tax_18S <-as.matrix(read.table("/PATH/TO/WORKING_DIRECTORY/ASVs_taxonomy.tsv", row.names = 1, header = TRUE, sep = "\t"))
map_18S <- read.table("/PATH/TO/WORKING_DIRECTORY/metadata.txt", row.names = 1, sep ="\t", header = TRUE) # This is your metadata file
ps_18S = phyloseq(otu_table(SV_18S, taxa_are_rows=TRUE), 
                  sample_data(map_18S), 
                  tax_table(tax_18S))
ps_18S
```
- You now have a phyloseq which is an R object containing all your data which you can filter and plot with
- Next step is to remove any non-protist contamination from your dataset. This includes any household plant and metazoan reads that may have been amplified during PCR.
- We will also remove any ASVs that have less than 200 reads total in our dataset, as these tend to be spurious and not important for our analysis.
- We will also remove reads that are not in the top 99% of at least 2 samples for the same reason. 
```{r}
psf_18S <- subset_taxa(ps_18S, (Subdivision!="Metazoa"| is.na(Division))) # This removes any metazoans. Hopefully not many since we used anti-metazoan primers
psf_18S <- subset_taxa(psf_18S, (Class!="Embryophyceae"| is.na(Class))) # This removes a common microbial contaminant. This is a land based family and should not be in our data.
ps_18S
psf_18S
psf_18S <- prune_samples(sample_sums(psf_18S)>=200, psf_18S)
psf_18S
f1<- filterfun_sample(topf(0.99))
wh1 <- genefilter_sample(psf_18S, f1, A=2)
psf_18S <- prune_taxa(wh1, psf_18S)
psf_18S
```

- This is to handle the NAs in our taxonomy and make sure we account for everything in our dataset, even those sequences not well annotated.
```{r}
library(stringr)
tax.clean <- psf_18S@tax_table
tax.clean <- data.frame(row.names = row.names(psf_18S@tax_table),
                        Domain = str_replace(psf_18S@tax_table[,1], "D_0__",""),
                        Supergroup = str_replace(psf_18S@tax_table[,2], "D_1__",""),
                        Division = str_replace(psf_18S@tax_table[,3], "D_2__",""),
                        Subdivision = str_replace(psf_18S@tax_table[,4], "D_3__",""),
                        Class = str_replace(psf_18S@tax_table[,5], "D_4__",""),
                        Order = str_replace(psf_18S@tax_table[,6], "D_5__",""),
                        Family = str_replace(psf_18S@tax_table[,7], "D_6__",""),
                        Genus = str_replace(psf_18S@tax_table[,8], "D_7__",""),
                        Species = str_replace(psf_18S@tax_table[,9], "D_8__",""),
                        stringsAsFactors = FALSE)
tax.clean[is.na(tax.clean)] <- ""

for (i in 1:9){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fille holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  
  if (tax.clean[i,2] == ""){domain <- paste(tax.clean[i,1],"_X", sep = "")
  tax.clean[i, 2] <- domain
  tax.clean[i, 3] <- paste(domain,"X", sep = "")
  tax.clean[i, 4] <- paste(domain,"XX", sep = "")
  tax.clean[i, 5] <- paste(domain,"XXX", sep = "")
  tax.clean[i, 6] <- paste(domain,"XXXX", sep = "")
  tax.clean[i, 7] <- paste(domain,"XXXXX", sep = "")
  tax.clean[i, 8] <- paste(domain,"XXXXXX", sep = "")
  tax.clean[i, 9] <- paste(domain,"XXXXXX_sp.", sep = "")
  } else 
    if (tax.clean[i,3] == ""){supergroup <- paste(tax.clean[i,2], "_X", sep = "")
    tax.clean[i, 3] <- supergroup
    tax.clean[i, 4] <- paste(supergroup,"X", sep = "")
    tax.clean[i, 5] <- paste(supergroup,"XX", sep = "")
    tax.clean[i, 6] <- paste(supergroup,"XXX", sep = "")
    tax.clean[i, 7] <- paste(supergroup,"XXXX.", sep = "")
    tax.clean[i, 8] <- paste(supergroup,"XXXXX", sep = "")
    tax.clean[i, 9] <- paste(supergroup,"XXXXX_sp.", sep = "")
    } else 
      if (tax.clean[i,4] == ""){division <- paste(tax.clean[i,3], "_X", sep = "")
      tax.clean[i, 4] <- division
      tax.clean[i, 5] <- paste(division,"X", sep = "")
      tax.clean[i, 6] <- paste(division,"XX", sep = "")
      tax.clean[i, 7] <- paste(division,"XXX", sep = "")
      tax.clean[i, 8] <- paste(division,"XXXX", sep = "")
      tax.clean[i, 9] <- paste(division,"XXXX_sp.", sep = "")
      } else 
        if (tax.clean[i,5] == ""){subdivision <- paste(tax.clean[i,4], "_X", sep = "")
        tax.clean[i, 5] <- subdivision
        tax.clean[i, 6] <- paste(subdivision,"X", sep = "")
        tax.clean[i, 7] <- paste(subdivision,"XX", sep = "")
        tax.clean[i, 8] <- paste(subdivision,"XXX", sep = "")
        tax.clean[i, 9] <- paste(subdivision,"XXX_sp.", sep = "")
        } else 
          if (tax.clean[i,6] == ""){class <- paste(tax.clean[i,5], "_X", sep = "")
          tax.clean[i, 6] <- class
          tax.clean[i, 7] <- paste(class,"X", sep = "")
          tax.clean[i, 8] <- paste(class,"XX", sep = "")
          tax.clean[i, 9] <- paste(class,"XX_sp.", sep = "")
          } else 
            if (tax.clean[i,7] == ""){order <- paste(tax.clean[i,6], "_X", sep = "")
            tax.clean[i, 7] <- order
            tax.clean[i, 8] <- paste(order,"X", sep = "")
            tax.clean[i, 9] <- paste(order,"X_sp.", sep = "")
            } else 
              if (tax.clean[i,8] == ""){family <- paste(tax.clean[i,7], "_X", sep = "")
              tax.clean[i, 8] <- family
              tax.clean[i, 9] <- paste(family,"X_sp.", sep = "")
              } else 
                if (tax.clean[i,9] == ""){tax.clean$Species[i] <- paste(tax.clean$Genus[i], "sp.",sep = "_")}
}
View(tax.clean)
tax_table(psf_18S) <- as.matrix(tax.clean)
```
- now save your filtered phyloseq object for easy access next time. 
```{r}
saveRDS(psf_18S, file = "/PATH/TO/WORKING_DIRECTORY/psf_18S.rds")
```

- *** You can now load in the object and start here in the future! ***
  ```{r}
psf_18S <- readRDS("/Users/ellawilson/Desktop/Results/psf_18S.rds")
```

- Those are the basics to get you started! 
  - Now see what figures you can make using the code in the pdf!
  ```{r}
library(ALDEx2);packageVersion("ALDEx2")
library(vegan); packageVersion("vegan")
library(CoDaSeq); packageVersion("CoDaSeq")
taxon <- psf_18S@tax_table
d.czm <- cmultRepl(psf_18S@otu_table, method="CZM", label=0)
d.clr <- codaSeq.clr(d.czm)
E.clr <- t(d.clr)
d.pcx <- prcomp(E.clr)
dist.clr <- dist(E.clr)
var1 <- psf_18S@sam_data$Heat.Stress
ano.var1 <- anosim(dist.clr, var1, permutations=999)
plot(ano.var1)

```
```{r}
library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)
### Adjust below code "SD" and what not accordingly.
### VAR1 is the variable you are interested in. VAR2 (or VAR3, VAR4, etcs) are
output_combined = ancombc2(data = psf_18S, assay_name = "counts", tax_level =
                             "Genus",fix_formula = "Heat.Stress", rand_formula = NULL,p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,group = "Heat.Stress", struc_zero = TRUE, neg_lb = FALSE,alpha = 0.05, n_cl = 2, verbose = TRUE,global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))
res_prim = output_combined$res_pair
res_prim
write.csv(res_prim,"/Users/ellawilson/Desktop/Results/ANCOM.CSV")
```


```{r}
library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)
### Adjust below code "SD" and what not accordingly.
### VAR1 is the variable you are interested in. VAR2 (or VAR3, VAR4, etcs) are
output_combined = ancombc2(data = psf_18S, assay_name = "counts", tax_level =
                             "Family",fix_formula = "Heat.Stress", rand_formula = NULL,p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,group = "Heat.Stress", struc_zero = TRUE, neg_lb = FALSE,alpha = 0.05, n_cl = 2, verbose = TRUE,global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))
res_prim = output_combined$res_pair
res_prim
write.csv(res_prim,"/Users/ellawilson/Desktop/Results/ANCOM_FAMILY.CSV")
```


```{r}
library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)
### Adjust below code "SD" and what not accordingly.
### VAR1 is the variable you are interested in. VAR2 (or VAR3, VAR4, etcs) are
output_combined = ancombc2(data = psf_18S, assay_name = "counts", tax_level =
                             NULL,fix_formula = "Heat.Stress", rand_formula = NULL,p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,group = "Heat.Stress", struc_zero = TRUE, neg_lb = FALSE,alpha = 0.05, n_cl = 2, verbose = TRUE,global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))
res_prim = output_combined$res_pair
res_prim
write.csv(res_prim,"/Users/ellawilson/Desktop/Results/ANCOM_OTU.CSV")
```


```{r}
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)
psf_18S_nocold
psf_18S
sample_data(psf_18S_nocold)$Heat.Stress
sample_data(psf_18S_nocold)$Heat.Stress <- factor(
  sample_data(psf_18S_nocold)$Heat.Stress,
  levels = c("control", "mid way", "heat"), #Or whatever the variables are called in that column
  ordered = TRUE
)

library(ANCOMBC)
library(tidyverse)

output_combined <- ancombc2(
  data = psf_18S_nocold,
  assay_name = "counts",
  tax_level = "Family",
  
  fix_formula = "Heat.Stress",
  rand_formula = NULL,
  
  p_adj_method = "fdr",
  
  pseudo = 0,
  pseudo_sens = TRUE,
  
  prv_cut = 0.01,
  lib_cut = 0,
  s0_perc = 0.05,
  
  group = "Heat.Stress",
  
  struc_zero = TRUE,
  neg_lb = FALSE,
  
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  
  ## 🔥 TREND TEST
  trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_FAMILY_trend.csv"
)
```

```{r}
psf_18S <- readRDS("/Users/ellawilson/Desktop/Results/psf_18S.rds")
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)

psf_18S_nocold@sam_data$Heat.Stress <- factor(
  psf_18S_nocold@sam_data$Heat.Stress,
  levels = c("control", "mid way", "heat"),
  labels = c("A-control", "B-midway", "C-heat")
)

library(ANCOMBC)
library(tidyverse)
library(DT)
output_combined <- ancombc2(data = psf_18S_nocold,
                            assay_name = "counts",
                            tax_level = "Family",
                            fix_formula = "Heat.Stress",
                            rand_formula = NULL,
                            p_adj_method = "fdr",
                            pseudo = 0,
                            pseudo_sens = TRUE,
                            prv_cut = 0.01,
                            lib_cut = 0,
                            s0_perc = 0.05,
                            group = "Heat.Stress",
                            struc_zero = TRUE,
                            neg_lb = FALSE,
                            alpha = 0.05,
                            n_cl = 2,
                            verbose = TRUE,
                            global = TRUE,
                            pairwise = FALSE,
                            dunnet = FALSE,
                            trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE)),
                                                 node = list(2, 2, 1),
                                                 solver = "ECOS",
                                                 B = 10))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_Family_trend.csv"
)
```
```{r}
psf_18S <- readRDS("/Users/ellawilson/Desktop/Results/psf_18S.rds")
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)

psf_18S_nocold@sam_data$Heat.Stress <- factor(
  psf_18S_nocold@sam_data$Heat.Stress,
  levels = c("control", "mid way", "heat"),
  labels = c("A-control", "B-midway", "C-heat")
)

library(ANCOMBC)
library(tidyverse)
library(DT)
output_combined <- ancombc2(data = psf_18S_nocold,
                            assay_name = "counts",
                            tax_level = "Genus",
                            fix_formula = "Heat.Stress",
                            rand_formula = NULL,
                            p_adj_method = "fdr",
                            pseudo = 0,
                            pseudo_sens = TRUE,
                            prv_cut = 0.01,
                            lib_cut = 0,
                            s0_perc = 0.05,
                            group = "Heat.Stress",
                            struc_zero = TRUE,
                            neg_lb = FALSE,
                            alpha = 0.05,
                            n_cl = 2,
                            verbose = TRUE,
                            global = TRUE,
                            pairwise = FALSE,
                            dunnet = FALSE,
                            trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE)),
                                                 node = list(2, 2, 1),
                                                 solver = "ECOS",
                                                 B = 10))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_Genus_trend.csv"
)
```
```{r}
psf_18S <- readRDS("/Users/ellawilson/Desktop/Results/psf_18S.rds")
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)

psf_18S_nocold@sam_data$Heat.Stress <- factor(
  psf_18S_nocold@sam_data$Heat.Stress,
  levels = c("control", "mid way", "heat"),
  labels = c("A-control", "B-midway", "C-heat")
)

library(ANCOMBC)
library(tidyverse)
library(DT)
output_combined <- ancombc2(data = psf_18S_nocold,
                            assay_name = "counts",
                            tax_level = "Species",
                            fix_formula = "Heat.Stress",
                            rand_formula = NULL,
                            p_adj_method = "fdr",
                            pseudo = 0,
                            pseudo_sens = TRUE,
                            prv_cut = 0.01,
                            lib_cut = 0,
                            s0_perc = 0.05,
                            group = "Heat.Stress",
                            struc_zero = TRUE,
                            neg_lb = FALSE,
                            alpha = 0.05,
                            n_cl = 2,
                            verbose = TRUE,
                            global = TRUE,
                            pairwise = FALSE,
                            dunnet = FALSE,
                            trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE)),
                                                 node = list(2, 2, 1),
                                                 solver = "ECOS",
                                                 B = 10))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_Species_trend.csv"
)
```
```{r}
psf_18S <- readRDS("/Users/ellawilson/Desktop/Results/psf_18S.rds")
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)

psf_18S_nocold@sam_data$Heat.Stress <- factor(
  psf_18S_nocold@sam_data$Heat.Stress,
  levels = c("control", "mid way", "heat"),
  labels = c("A-control", "B-midway", "C-heat")
)

library(ANCOMBC)
library(tidyverse)
library(DT)
output_combined <- ancombc2(data = psf_18S_nocold,
                            assay_name = "counts",
                            tax_level = NULL,
                            fix_formula = "Heat.Stress",
                            rand_formula = NULL,
                            p_adj_method = "fdr",
                            pseudo = 0,
                            pseudo_sens = TRUE,
                            prv_cut = 0.01,
                            lib_cut = 0,
                            s0_perc = 0.05,
                            group = "Heat.Stress",
                            struc_zero = TRUE,
                            neg_lb = FALSE,
                            alpha = 0.05,
                            n_cl = 2,
                            verbose = TRUE,
                            global = TRUE,
                            pairwise = FALSE,
                            dunnet = FALSE,
                            trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE)),
                                                 node = list(2, 2, 1),
                                                 solver = "ECOS",
                                                 B = 10))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_NULL_trend.csv"
)
```
```{r}
psf_18S <- readRDS("/Users/ellawilson/Desktop/Results/psf_18S.rds")
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)

psf_18S_nocold@sam_data$Heat.Stress <- factor(
  psf_18S_nocold@sam_data$Heat.Stress,
  levels = c("control", "mid way", "heat"),
  labels = c("A-control", "B-midway", "C-heat")
)

library(ANCOMBC)
library(tidyverse)
library(DT)
output_combined <- ancombc2(data = psf_18S_nocold,
                            assay_name = "counts",
                            tax_level = "Genus",
                            fix_formula = "Heat.Stress",
                            rand_formula = NULL,
                            p_adj_method = "fdr",
                            pseudo = 0,
                            pseudo_sens = TRUE,
                            prv_cut = 0.01,
                            lib_cut = 0,
                            s0_perc = 0.05,
                            group = "Heat.Stress",
                            struc_zero = TRUE,
                            neg_lb = FALSE,
                            alpha = 0.05,
                            n_cl = 2,
                            verbose = TRUE,
                            global = TRUE,
                            pairwise = FALSE,
                            dunnet = FALSE,
                            trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),
                            em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),
                            mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
                            trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(-1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE),
                                                                 matrix(c(1, 0, 1, -1),
                                                                        nrow = 2, 
                                                                        byrow = TRUE)),
                                                 node = list(2, 2, 1),
                                                 solver = "ECOS",
                                                 B = 10))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_Genus_trend.csv"
)

saveRDS(res_trend, "/Users/ellawilson/Desktop/Results/res_trend.rds")
res_trend <- readRDS("/Users/ellawilson/Desktop/Results/res_trend.rds")
library(dplyr)

df_fig_trend <- res_trend %>%
  dplyr::filter(diff_abn == 1) %>%
  dplyr::mutate(
    lfc1 = round(`lfc_Heat.StressB-midway`, 2),
    lfc2 = round(`lfc_Heat.StressC-heat`, 2),
    color = ifelse(diff_abn, "aquamarine3", "black")
  ) %>%
  tidyr::pivot_longer(
    cols = c(lfc1, lfc2),
    names_to = "group",
    values_to = "value"
  ) %>%
  dplyr::arrange(taxon)


df_fig_trend$group = recode(df_fig_trend$group, 
                            lfc1 = "Midway - Control",
                            lfc2 = "Heat - Control")
df_fig_trend$group = factor(df_fig_trend$group, 
                            levels = c("Midway - Control", 
                                       "Heat - Control"))

lo = floor(min(df_fig_trend$value))
up = ceiling(max(df_fig_trend$value))
mid = (lo + up)/2
library(ggplot2)
fig_trend = df_fig_trend %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to control coral") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color = df_fig_trend %>%
                                     dplyr::distinct(taxon, color) %>%
                                     .$color))
fig_trend
```

```{r}
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)
psf_18S_nocold
psf_18S
sample_data(psf_18S_nocold)$Heat.Stress
sample_data(psf_18S_nocold)$Heat.Stress <- factor(
  sample_data(psf_18S_nocold)$Heat.Stress,
  levels = c("control", "mid way", "heat"), #Or whatever the variables are called in that column
  ordered = TRUE
)

library(ANCOMBC)
library(tidyverse)

output_combined <- ancombc2(
  data = psf_18S_nocold,
  assay_name = "counts",
  tax_level = "Genus",
  
  fix_formula = "Heat.Stress",
  rand_formula = NULL,
  
  p_adj_method = "fdr",
  
  pseudo = 0,
  pseudo_sens = TRUE,
  
  prv_cut = 0.01,
  lib_cut = 0,
  s0_perc = 0.05,
  
  group = "Heat.Stress",
  
  struc_zero = TRUE,
  neg_lb = FALSE,
  
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  
  ## 🔥 TREND TEST
  trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_GENUS_trend.csv"
)
```
```{r}
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)
psf_18S_nocold
psf_18S
sample_data(psf_18S_nocold)$Heat.Stress
sample_data(psf_18S_nocold)$Heat.Stress <- factor(
  sample_data(psf_18S_nocold)$Heat.Stress,
  levels = c("control", "mid way", "heat"), #Or whatever the variables are called in that column
  ordered = TRUE
)

library(ANCOMBC)
library(tidyverse)

output_combined <- ancombc2(
  data = psf_18S_nocold,
  assay_name = "counts",
  tax_level = "Species",
  
  fix_formula = "Heat.Stress",
  rand_formula = NULL,
  
  p_adj_method = "fdr",
  
  pseudo = 0,
  pseudo_sens = TRUE,
  
  prv_cut = 0.01,
  lib_cut = 0,
  s0_perc = 0.05,
  
  group = "Heat.Stress",
  
  struc_zero = TRUE,
  neg_lb = FALSE,
  
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  
  ## 🔥 TREND TEST
  trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_SPECIES_trend.csv"
)
```
```{r}
psf_18S_nocold <- subset_samples(
  psf_18S,
  Heat.Stress != "cold"
)
psf_18S_nocold
psf_18S
sample_data(psf_18S_nocold)$Heat.Stress
sample_data(psf_18S_nocold)$Heat.Stress <- factor(
  sample_data(psf_18S_nocold)$Heat.Stress,
  levels = c("control", "mid way", "heat"), #Or whatever the variables are called in that column
  ordered = TRUE
)

library(ANCOMBC)
library(tidyverse)

output_combined <- ancombc2(
  data = psf_18S_nocold,
  assay_name = "counts",
  tax_level = NULL,
  
  fix_formula = "Heat.Stress",
  rand_formula = NULL,
  
  p_adj_method = "fdr",
  
  pseudo = 0,
  pseudo_sens = TRUE,
  
  prv_cut = 0.01,
  lib_cut = 0,
  s0_perc = 0.05,
  
  group = "Heat.Stress",
  
  struc_zero = TRUE,
  neg_lb = FALSE,
  
  alpha = 0.05,
  n_cl = 2,
  verbose = TRUE,
  
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  
  ## 🔥 TREND TEST
  trend = TRUE,iter_control = list(tol = 1e-2, max_iter = 20,verbose = TRUE),em_control = list(tol = 1e-5, max_iter = 100),lme_control = lme4::lmerControl(),mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),nrow = 2,byrow = TRUE),matrix(c(-1, 0, 1, -1),nrow = 2,byrow = TRUE)),node = list(2, 2),solver = "ECOS",B = 100))

res_trend <- output_combined$res_trend
res_trend
write.csv(
  res_trend,
  "/Users/ellawilson/Desktop/Results/ANCOM_NULL_trend.csv"
)
```


```

```{r}
res_prim = output_combined$res_pair
res_prim
write.csv(res_prim,"/Users/ellawilson/Desktop/Results/ANCOM.CSV")
df_res = res_prim %>% dplyr::filter(diff_Heat.Stressheat_Heat.Stresscontrol) %>%
  dplyr::select(taxon, contains("Heat.Stress"))
df_fig_res = df_res %>%arrange(desc(lfc_Heat.Stressheat_Heat.Stresscontrol)) %>%mutate(direct = ifelse(lfc_Heat.Stressheat_Heat.Stresscontrol> 0, "Positive LFC", "Negative LFC"))
df_fig_res$taxon = factor(df_fig_res$taxon, levels = df_fig_res$taxon)
df_fig_res$direct = factor(df_fig_res$direct,levels = c("Positive LFC", "Negative LFC"))
fig_res = df_fig_res %>%ggplot(aes(x = taxon, y = lfc_Heat.Stressheat_Heat.Stresscontrol, fill = lfc_Heat.Stressheat_Heat.Stresscontrol)) + geom_bar(stat = "identity", width = 0.7, color = "black",position = position_dodge(width = 0.4)) + geom_errorbar(aes(ymin = lfc_Heat.Stressheat_Heat.Stresscontrol - se_Heat.Stressheat_Heat.Stresscontrol, ymax = lfc_Heat.Stressheat_Heat.Stresscontrol +se_Heat.Stressheat_Heat.Stresscontrol),width = 0.2, position = position_dodge(0.05), color = "black") +labs(x = "Prokaryotic Order", y = "Log fold change",title = "LFC in Heat vs Control Corals", subtitle = "Of ANCOM-BC signficant taxa") + scale_fill_viridis_c(option = "plasma") +scale_color_viridis_c(option = "plasma") +theme_bw() +theme(legend.position = "none", panel.grid.minor.y = element_blank()) + coord_flip()
fig_res
```

```{r}
library(phyloseq)
library(dplyr)
library(ggplot2)
library(data.table)
```

````{r}
tax.sort <- as.data.frame(psf_18S@tax_table)
tax.sort <- tax.sort %>%arrange(Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species)
euk_genera <- tax.sort %>% pull(Genus) %>%unique() %>%toString()
# Split the string into individual values
euk_genera <- strsplit(euk_genera, ", ")[[1]]
euk_genera <- c(euk_genera, "Below Threshold")
sd1 <- psf_18S %>%tax_glom(taxrank = "Genus") %>%transform_sample_counts(function(x) { x / sum(x) * 100 }) %>%psmelt() %>%group_by(OTU, Heat.Stress) %>%summarize(sd_within = sd(Abundance, na.rm = FALSE),.groups = "drop") # keep structure for join
bub1 <- merge_samples(psf_18S, "Heat.Stress")
type_taxa1 <- bub1 %>% tax_glom(taxrank ="Genus") %>% transform_sample_counts(function(x) {x/sum(x)*100} )
type_taxa1 <- type_taxa1 %>% psmelt() %>% mutate(Heat.Stress = as.character(Sample)) %>% arrange(OTU)
library(data.table)
# create dataframe from phyloseq object
dat1 <- data.table(type_taxa1)  %>%left_join(sd1, by = c("OTU", "Heat.Stress")) %>%as.data.table()
# convert Phylum to a character vector from a factor because R
dat1$Genus <- as.character(dat1$Genus)
# group dataframe by Genus, calculate mean rel. abundance
dat1[, mean := mean(Abundance, na.rm = FALSE),by = "Genus"]
# Change name to remainder of Family less than 0.01%
dat1[(mean <= 0.01), Genus := "Below Threshold"]
dat1[(Genus == "Below Threshold"), Family := "Other"]
my_title1 <- expression(paste("Microbial Genera within ", italic("Pocillopora")))
Fig1A <- ggplot(dat1[Abundance > 0], aes(x = Sample, y = factor(Genus, levels=euk_genera))) + geom_point(aes(size = Abundance,fill= Family, color=Family), alpha = 0.5, shape = 21) +scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks =c(1,10,50,75)) + labs( x= "Thermal Stress", y = "Genus", size = "Relative Abundance %", fill="Family") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +ggtitle(label=my_title1, subtitle = "With mean relative abundance above 0.01%") + geom_point(aes(size=Abundance+sd_within, fill = Family, color=Family),alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-sd_within, fill = Family,color=Family), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
Fig1A

```
```{r}
c25 <- c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00","black", "gold1","skyblue2", "#FB9A99", "palegreen2","#CAB2D6","#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1","steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

psf_18S_nozoox <- subset_taxa(psf_18S, (Family!="Symbiodiniaceae"| is.na(Family))) # This removes any metazoans. Hopefully not many since we used anti-metazoan primers```
tax.sort <- as.data.frame(psf_18S_nozoox@tax_table)
tax.sort <- tax.sort %>%arrange(Domain, Supergroup, Division, Subdivision, Class, Order, Family, Genus, Species)
euk_genera <- tax.sort %>% pull(Genus) %>%unique() %>%toString()
# Split the string into individual values
euk_genera <- strsplit(euk_genera, ", ")[[1]]
euk_genera <- c(euk_genera, "Below Threshold")
samples_ordered <- c("cold", "control", "mid way", "heat")
sd1 <- psf_18S_nozoox %>%tax_glom(taxrank = "Genus") %>%transform_sample_counts(function(x) { x / sum(x) * 100 }) %>%psmelt() %>%group_by(OTU, Heat.Stress) %>%summarize(sd_within = sd(Abundance, na.rm = FALSE),.groups = "drop") # keep structure for join
bub1 <- merge_samples(psf_18S_nozoox, "Heat.Stress")
type_taxa1 <- bub1 %>% tax_glom(taxrank ="Genus") %>% transform_sample_counts(function(x) {x/sum(x)*100} )
type_taxa1 <- type_taxa1 %>% psmelt() %>% mutate(Heat.Stress = as.character(Sample)) %>% arrange(OTU)
library(data.table)
# create dataframe from phyloseq object
dat1 <- data.table(type_taxa1)  %>%left_join(sd1, by = c("OTU", "Heat.Stress")) %>%as.data.table()
# convert Phylum to a character vector from a factor because R
dat1$Genus <- as.character(dat1$Genus)
# group dataframe by Genus, calculate mean rel. abundance
dat1[, mean := mean(Abundance, na.rm = FALSE),by = "Genus"]
# Change name to remainder of Family less than 0.01%
dat1[(mean <= 0.01), Genus := "Below Threshold"]
dat1[(Genus == "Below Threshold"), Subdivision := "Other"]
my_title1 <- expression(paste("Microbial Genera within ", italic("Pocillopora")))
Fig1A <- ggplot(dat1[Abundance > 0], aes(x = factor(Sample, levels= samples_ordered), y = factor(Genus, levels=euk_genera))) + geom_point(aes(size = Abundance,fill= Subdivision, color=Subdivision), alpha = 0.5, shape = 21) +scale_size_continuous(limits = c(0.000001, 100), range = c(0.1,10), breaks =c(1,10,50,75)) + labs( x= "Thermal Stress", y = "Genus", size = "Relative Abundance %", fill="Subdivision") + theme_bw() + scale_fill_manual(values = c25) + scale_color_manual(values = c25) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust=1)) +ggtitle(label=my_title1, subtitle = "With mean relative abundance above 0.01%") + geom_point(aes(size=Abundance+sd_within, fill = Subdivision, color=Subdivision),alpha = 0.20, shape=21) + geom_point(aes(size=Abundance-sd_within, fill = Subdivision,color=Subdivision), alpha = 1, shape=21) + guides(color = FALSE, fill = guide_legend(override.aes = list(size = 3)))
Fig1A

```
```{r}
ggsave(filename ="/Users/ellawilson/Desktop/Results/bubble_plot.png",dpi=300)