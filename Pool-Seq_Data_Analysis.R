# R Analysis Script for Chironomus riparius Pool-Seq Data
# Author: Cosima Caliendof
# Date: January 2021

# Load required packages
library(dplyr)

#################################################
# FET and FST Analysis Using Popoolation2 Output
#################################################

# Set working directory
setwd("/path/to/analysis/directory")

# Read in FET table (contains log p-values in columns $4-$9)
rot_allGen_fet <- read.table("chiro_adaptation_pool_rot_allGen.fet")

# Calculate actual p-values (10^-logp) for each value
rot_allGen_fet$qval1 <- 10^-rot_allGen_fet$V4
rot_allGen_fet$qval2 <- 10^-rot_allGen_fet$V5
rot_allGen_fet$qval3 <- 10^-rot_allGen_fet$V6
rot_allGen_fet$qval4 <- 10^-rot_allGen_fet$V7
rot_allGen_fet$qval5 <- 10^-rot_allGen_fet$V8
rot_allGen_fet$qval6 <- 10^-rot_allGen_fet$V9

# Apply Benjamini-Hochberg correction for multiple testing
rot_allGen_fet$BHcorr1 <- p.adjust(rot_allGen_fet$qval1, method = "fdr")
rot_allGen_fet$BHcorr2 <- p.adjust(rot_allGen_fet$qval2, method = "fdr")
rot_allGen_fet$BHcorr3 <- p.adjust(rot_allGen_fet$qval3, method = "fdr")
rot_allGen_fet$BHcorr4 <- p.adjust(rot_allGen_fet$qval4, method = "fdr")
rot_allGen_fet$BHcorr5 <- p.adjust(rot_allGen_fet$qval5, method = "fdr")
rot_allGen_fet$BHcorr6 <- p.adjust(rot_allGen_fet$qval6, method = "fdr")

# Read in FST table (contains FST values in columns $3-$8)
rot_allGen_fst <- read.table("chiro_adaptation_pool_rot_allGen.fst")

# Extract relevant columns for analysis
rot_allGen_fst_1u2 <- rot_allGen_fst[,1:3]
rot_allGen_fet_1u2 <- subset(rot_allGen_fet, select = c(V1, V2, V4, qval1, BHcorr1))

# Filter for FST values > 0.1
rot_allGen_fst_1u2_filtered <- rot_allGen_fst[rot_allGen_fst$V3 >= 0.1, 1:3]

# Merge FET and FST results
# 1. Complete merge
rot_allGen_FET_FST_1u2 <- merge(rot_allGen_fst_1u2, rot_allGen_fet_1u2, by = c("V1","V2"), all = FALSE)
# 2. Filter for significant FET values (BH-corrected p < 0.01)
rot_allGen_FET_FST_1u2_BH <- subset(rot_allGen_FET_FST_1u2, BHcorr1 <= 0.01)

# 3. Merge filtered FST with FET
rot_allGen_FET_FST_1u2_filtered <- merge(rot_allGen_fst_1u2_filtered, rot_allGen_fet_1u2, by = c("V1","V2"), all = FALSE)
# 4. Filter for significant FET values (BH-corrected p < 0.01)
rot_allGen_FET_FST_1u2_filtered_BH <- subset(rot_allGen_FET_FST_1u2_filtered, BHcorr1 <= 0.01)

# Save results to files
write.csv(rot_allGen_FET_FST_1u2_BH, file = "rot_allGen_FET_FST_1u2_BH.csv")
write.csv(rot_allGen_FET_FST_1u2_filtered_BH, file = "rot_allGen_FET_FST_1u2_filtered_BH.csv")

#################################################
# Analysis with the poolSeq R Package
#################################################

# Load poolSeq package
library(poolSeq)

# Read sync file and extract treatment replicates
Pool_Rot <- read.sync(file="chiro_pool_rot.idf.sync", 
                     gen=c(1,2,3,4,5,6,7), 
                     repl=c(1,1,1,1,1,1,1),
                     keepOnlyBiallelic=TRUE)

# Extract allele frequencies for each generation
af_rot_gen1 <- af(sync=Pool_Rot, repl=1, gen=1)
af_rot_gen2 <- af(sync=Pool_Rot, repl=1, gen=2)
af_rot_gen3 <- af(sync=Pool_Rot, repl=1, gen=3)
af_rot_gen4 <- af(sync=Pool_Rot, repl=1, gen=4)
af_rot_gen5 <- af(sync=Pool_Rot, repl=1, gen=5)
af_rot_gen6 <- af(sync=Pool_Rot, repl=1, gen=6)
af_rot_gen7 <- af(sync=Pool_Rot, repl=1, gen=7)

# Save allele frequencies to files
write.table(af_rot_gen1, "af_rot_gen1.txt", sep="\t")
write.table(af_rot_gen2, "af_rot_gen2.txt", sep="\t")
write.table(af_rot_gen3, "af_rot_gen3.txt", sep="\t")
write.table(af_rot_gen4, "af_rot_gen4.txt", sep="\t")
write.table(af_rot_gen5, "af_rot_gen5.txt", sep="\t")
write.table(af_rot_gen6, "af_rot_gen6.txt", sep="\t")
write.table(af_rot_gen7, "af_rot_gen7.txt", sep="\t")

# Calculate coverage for each generation
cov_rot_gen1 <- coverage(sync=Pool_Rot, repl=1, gen=1)
cov_rot_gen2 <- coverage(sync=Pool_Rot, repl=1, gen=2)
cov_rot_gen3 <- coverage(sync=Pool_Rot, repl=1, gen=3)
cov_rot_gen4 <- coverage(sync=Pool_Rot, repl=1, gen=4)
cov_rot_gen5 <- coverage(sync=Pool_Rot, repl=1, gen=5)
cov_rot_gen6 <- coverage(sync=Pool_Rot, repl=1, gen=6)
cov_rot_gen7 <- coverage(sync=Pool_Rot, repl=1, gen=7)

# Calculate allele frequency changes (AFC) between consecutive generations
AFC_rot_1vs2 <- af_rot_gen1 - af_rot_gen2
AFC_rot_2vs3 <- af_rot_gen2 - af_rot_gen3
AFC_rot_3vs4 <- af_rot_gen3 - af_rot_gen4
AFC_rot_4vs5 <- af_rot_gen4 - af_rot_gen5
AFC_rot_5vs6 <- af_rot_gen5 - af_rot_gen6
AFC_rot_6vs7 <- af_rot_gen6 - af_rot_gen7
AFC_rot_1vs7 <- af_rot_gen1 - af_rot_gen7  # Overall change from gen 1 to 7

# Save AFC values to files
write.table(AFC_rot_1vs2, "afc_rot_gen1vs2.txt", sep="\t")
write.table(AFC_rot_2vs3, "afc_rot_gen2vs3.txt", sep="\t")
write.table(AFC_rot_3vs4, "afc_rot_gen3vs4.txt", sep="\t")
write.table(AFC_rot_4vs5, "afc_rot_gen4vs5.txt", sep="\t")
write.table(AFC_rot_5vs6, "afc_rot_gen5vs6.txt", sep="\t")
write.table(AFC_rot_6vs7, "afc_rot_gen6vs7.txt", sep="\t")
write.table(AFC_rot_1vs7, "afc_rot_gen1vs7.txt", sep="\t")

# Prepare for Chi-square tests
# Create count matrices for each generation
# Calculate reference allele counts (number of reads supporting the reference allele)
NA_rot_gen1 <- t(cov_rot_gen1 * af_rot_gen1)
NA_rot_gen2 <- t(cov_rot_gen2 * af_rot_gen2)
NA_rot_gen3 <- t(cov_rot_gen3 * af_rot_gen3)
NA_rot_gen4 <- t(cov_rot_gen4 * af_rot_gen4)
NA_rot_gen5 <- t(cov_rot_gen5 * af_rot_gen5)
NA_rot_gen6 <- t(cov_rot_gen6 * af_rot_gen6)
NA_rot_gen7 <- t(cov_rot_gen7 * af_rot_gen7)

# Calculate alternative allele counts
na_rot_gen1 <- cov_rot_gen1 - NA_rot_gen1
na_rot_gen2 <- cov_rot_gen2 - NA_rot_gen2
na_rot_gen3 <- cov_rot_gen3 - NA_rot_gen3
na_rot_gen4 <- cov_rot_gen4 - NA_rot_gen4
na_rot_gen5 <- cov_rot_gen5 - NA_rot_gen5
na_rot_gen6 <- cov_rot_gen6 - NA_rot_gen6
na_rot_gen7 <- cov_rot_gen7 - NA_rot_gen7

# Run Chi-square tests for consecutive generations
# Parameter explanation:
#   min.cov: Minimum coverage threshold
#   min.cnt: Minimum count threshold
#   max.cov: Maximum coverage threshold
#   log: Return -log10(p-values) if TRUE
p.values_delta1 <- chi.sq.test(A0=NA_rot_gen1, a0=na_rot_gen1, At=NA_rot_gen2, at=na_rot_gen2, 
                              min.cov=15, min.cnt=3, max.cov=120, log=TRUE)
p.values_delta2 <- chi.sq.test(A0=NA_rot_gen2, a0=na_rot_gen2, At=NA_rot_gen3, at=na_rot_gen3, 
                              min.cov=15, min.cnt=3, max.cov=120, log=TRUE)
p.values_delta3 <- chi.sq.test(A0=NA_rot_gen3, a0=na_rot_gen3, At=NA_rot_gen4, at=na_rot_gen4, 
                              min.cov=15, min.cnt=3, max.cov=120, log=TRUE)
p.values_delta4 <- chi.sq.test(A0=NA_rot_gen4, a0=na_rot_gen4, At=NA_rot_gen5, at=na_rot_gen5, 
                              min.cov=15, min.cnt=3, max.cov=120, log=TRUE)
p.values_delta5 <- chi.sq.test(A0=NA_rot_gen5, a0=na_rot_gen5, At=NA_rot_gen6, at=na_rot_gen6, 
                              min.cov=15, min.cnt=3, max.cov=120, log=TRUE)
p.values_delta6 <- chi.sq.test(A0=NA_rot_gen6, a0=na_rot_gen6, At=NA_rot_gen7, at=na_rot_gen7, 
                              min.cov=15, min.cnt=3, max.cov=120, log=TRUE)

# Chi-square test for overall change (gen 1 vs gen 7)
p.values_delta1vs7 <- chi.sq.test(A0=NA_rot_gen1, a0=na_rot_gen1, At=NA_rot_gen7, at=na_rot_gen7, 
                                 min.cov=15, min.cnt=3, max.cov=120, log=TRUE)

# Calculate significance thresholds for each test
# Output shows the -log10(p-value) thresholds at different significance levels
quantile(p.values_delta1, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)
quantile(p.values_delta2, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)
quantile(p.values_delta3, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)
quantile(p.values_delta4, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)
quantile(p.values_delta5, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)
quantile(p.values_delta6, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)
quantile(p.values_delta1vs7, probs=c(0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)

# Save Chi-square test results
write.table(p.values_delta1, "p.values_delta1.txt", sep="\t")
write.table(p.values_delta2, "p.values_delta2.txt", sep="\t")
write.table(p.values_delta3, "p.values_delta3.txt", sep="\t")
write.table(p.values_delta4, "p.values_delta4.txt", sep="\t")
write.table(p.values_delta5, "p.values_delta5.txt", sep="\t")
write.table(p.values_delta6, "p.values_delta6.txt", sep="\t")
write.table(p.values_delta1vs7, "p.values_rot_delta1vs7.txt", sep="\t")

# Create Manhattan plots for each test
for(i in 1:6) {
  nam <- paste0("p.values_delta", i)
  mypath <- file.path("plots", paste("chiplot_Rot_", i, ".pdf", sep=""))
  pdf(file=mypath)
  plot(get(nam), 
       main=paste0("Chi-squared test - Replicate Rot, Gen", i, " vs ", (i+1)), 
       ylim=c(0, 20), 
       xlab="Positions", 
       ylab="-log10(p)", 
       pch=".")
  dev.off()
}

# Create Manhattan plot for overall change (gen 1 vs gen 7)
pdf("chiplot_Rot_1vs7.pdf")
plot(p.values_delta1vs7, 
     main="Chi-squared test - Replicate Rot, Gen1 vs Gen7", 
     ylim=c(0, 20), 
     xlab="Position", 
     ylab="-log10(p)", 
     pch=".")
dev.off()

# Apply Benjamini-Hochberg correction to overall change results
# Note: Need to convert from log-transformed p-values first
p.values_delta1vs7_raw <- chi.sq.test(A0=NA_rot_gen1, a0=na_rot_gen1, At=NA_rot_gen7, at=na_rot_gen7, 
                                     min.cov=15, min.cnt=3, max.cov=120, log=FALSE)
p.values_delta1vs7_bh <- p.adjust(p.values_delta1vs7_raw, method="fdr")
p.values_delta1vs7_logbh <- -log10(p.values_delta1vs7_bh)

# Calculate thresholds for BH-corrected p-values
quantile(p.values_delta1vs7_logbh, probs=c(0.95, 0.99, 0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)

# Save BH-corrected p-values
write.table(p.values_delta1vs7_logbh, "p.values_rot_delta1vs7_bh.txt", sep="\t")

# Plot BH-corrected results
pdf('chiplot_p.valbh_Rot_1vs7.pdf')
plot(p.values_delta1vs7_logbh, 
     main="Chi-squared test (BH-corrected) - Replicate Rot Gen1vs7", 
     ylim=c(0, 15), 
     xlab="Position", 
     ylab="-log10(p)", 
     pch=".")
dev.off()

#################################################
# Genetic Drift Simulation
#################################################

# Estimate effective population size (Ne)
Ne_Rot <- estimateNe(p0=af_rot_gen1, 
                    pt=af_rot_gen7, 
                    cov0=cov_rot_gen1, 
                    covt=cov_rot_gen7, 
                    t=6,                # 6 generations of change
                    Ncensus=100,        # Census population size
                    poolSize=c(100, 100))  # Sequencing pool size

print(Ne_Rot)  # Estimated Ne value

# Simulate expected allele frequency trajectories under drift
exp_Rot <- wf.traj(p0=af_rot_gen1, 
                  Ne=as.numeric(Ne_Rot), 
                  t=c(1,2,3,4,5,6,7),  # 7 time points
                  s=0,                 # No selection
                  haploid=FALSE)

# Save expected allele frequencies
for(i in 1:7) {
  nam <- paste0("af_rot_exp_gen", i)
  assign(nam, exp_Rot[,i])
  write.table(get(nam), paste0(nam, ".txt"), sep="\t")
}

# Compare expected vs. observed frequencies using Chi-square test
# Generate count matrices for expected frequencies
NA_rot_exp_gen7 <- t(cov_rot_gen7 * exp_Rot[,7])
na_rot_exp_gen7 <- cov_rot_gen7 - NA_rot_exp_gen7
NA_rot_exp_gen1 <- t(cov_rot_gen1 * exp_Rot[,1])
na_rot_exp_gen1 <- cov_rot_gen1 - NA_rot_exp_gen1

# Chi-square test comparing observed gen1 to expected gen7
p.values_rot_exp_delta1vs7 <- chi.sq.test(A0=NA_rot_exp_gen1, 
                                         a0=na_rot_exp_gen1, 
                                         At=NA_rot_exp_gen7, 
                                         at=na_rot_exp_gen7, 
                                         min.cov=15, 
                                         min.cnt=3, 
                                         max.cov=120, 
                                         log=TRUE)

# Calculate thresholds for drift simulation significance
quantile(p.values_rot_exp_delta1vs7, probs=c(0.95, 0.99, 0.999, 0.9999, 0.99999), na.rm=TRUE, names=TRUE)

# Save drift simulation results
write.table(p.values_rot_exp_delta1vs7, "p.values_rot_exp_delta1vs7.txt", sep="\t")

# Plot drift simulation results
pdf("driftplot_Rot_exp_1vs7.pdf")
plot(p.values_delta1vs7, 
     main="Expected vs. Observed Allele Frequencies - Replicate Rot, Gen1vs7", 
     ylim=c(0, 20), 
     xlab="Position", 
     ylab="-log10(p)", 
     pch=".")
dev.off()
