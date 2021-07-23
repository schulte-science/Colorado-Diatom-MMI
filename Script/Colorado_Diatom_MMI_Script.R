#Project title: Colorado diatom MMI computation and analysis

#Authors: Nicholas O. Schulte and Sarah A. Spaulding

#Date: 22 July 2021

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#####Process diatom data prior to input into MMI

#--------------

#Extract metadata and diatom counts from Colorado WQCD Algal Database
#-Load database
db <- read.csv("Colorado WQCD Algal Database as of 021920.csv", header = TRUE, strip.white = TRUE)

#-Retain diatom columns, identified as all columns between 'BACILLARIOPHYTA' and 'CRYPTOPHYTA' (row 1 values)
db <- db[,c(1:7,(which(db[1,] == "BACILLARIOPHYTA")+1):(which(db[1,] == "CRYPTOPHYTA")-2))]

#--------------

#Update and re-format diatom taxonomy
#-For any diatom taxon with 'Alternative Nomenclature' (row 2 names), replace 'nomenclature' (row 1) with 'alternative nomenclature'
for(i in 8:ncol(db)){
  db[1,i] <- ifelse(db[2,i] == "", db[1,i], db[2,i])
}

#-Replace double spaces in all columns (in case of typos), all spaces in species names with underscores, and all instances of "sp." (and anything after)
db[,1:ncol(db)] <- lapply(db[,1:ncol(db)], gsub, pattern = "  ", replacement = "")
db[,8:ncol(db)] <- lapply(db[,8:ncol(db)], gsub, pattern = " ", replacement = "_")
db[,8:ncol(db)] <- lapply(db[,8:ncol(db)], gsub, pattern = "_sp\\..*", replacement = "")

#--------------

#Re-format and fix typos in count data
#-Replace blank counts with "0" and "," with ".".
db[,8:ncol(db)][db[,8:ncol(db)] == ""] <- 0
db[,8:ncol(db)] <- lapply(db[,8:ncol(db)], gsub, pattern = ",", replacement = ".")

#--------------

#Sum and remove counts from duplicate taxa
#-Temporarily extract diatom names and counts only (i.e., remove metadata) (allows for having columns with duplicate names, necessary for subsequent code)
db_tax <- as.matrix(as.data.frame(db[c(1,3:nrow(db)),8:ncol(db)]))
colnames(db_tax) <- db_tax[1,]
db_tax <- db_tax[-c(1),]

#-Retain CDPHE taxon names
cdphe_tax <- colnames(db_tax)

#-Update taxon names for consistency with taxonomy used in MMI Shiny application (https://diatom-mmi-demo.shinyapps.io/Diatom_MMI_V24/)
#--Load taxonomy lookup table manually developed as part of this project
#---Note that this lookup table considers CDPHE taxonomy from the WQCD Algal Database as of 01/2019 to the taxonomy used in the MMI as of 05/2021
lookup <- read.csv("CDPHE_Taxonomy_Lookup_Table.csv", header = TRUE, strip.white = TRUE)

#--Match CDPHE names in db_long to MMI names in lookup for db_long, db_rel, and db
overlap <- match(colnames(db_tax), lookup$CDPHE)
colnames(db_tax)[!is.na(overlap)] <- lookup$MMI[na.omit(overlap)]

rm(lookup, overlap)

#--Future revisions to the taxonomy lookup table may be required, based on changes to CDPHE and MMI taxonomy
#---Steps on how to check CDPHE taxonomy using the MMI Shiny application are outlined below
#----Upload "CDPHE_Sample_Taxonomy_Counts.csv" in the Shiny application - Step 2: Upload Taxonomic Data
#-----Be sure to click the dropdown menu in the bottom right to search for "All Files" instead of "Custom Files" to display the file
#----Download "MISSING TAXA LIST"
#----Compare taxonomy in "MISSING TAXA LIST" to BioData using 'biodata_check' (https://github.com/bishopia/biodata_check)
#-----'biodata_check' may not resolve taxonomic issues, as the taxonomy used by the MMI may vary, as BioData names may be outdated
#-----A future update to the MMI Shiny application will be a downloadable list of taxon names used in the MMI, against which you can check CDPHE names directly

#-Change matrix to numeric
db_tax <- apply(db_tax, 2 ,as.numeric)

#-Identify taxa with counts in multiple columns
dups <- colnames(db_tax[,which(duplicated(colnames(db_tax)))])
dups

#-Sum rows for each duplicated taxon name
for(a in dups){
  db_tax[,colnames(db_tax) == a] <- rowSums(db_tax[,colnames(db_tax) == a])
}

#-Remove duplicated taxon columns (each column with the same name has the same rowSum data)
db_tax <- db_tax[,-which(duplicated(colnames(db_tax)))]

#-Merge temporary diatom count matrix with main db
#--Create temporary dataframe for metadata and replace column names with second row
db_meta <- db[,1:7]
colnames(db_meta) <- db_meta[2,]
db_meta <- db_meta[-c(1:2),]

#--Merge db_meta and db_tax back into db, then remove all temporary objects
db <- cbind(db_meta, db_tax)
rm(db_meta, db_tax, a, dups, i)

#--------------

#Process metadata
#-Add unique sample identifier (UID) to each sample and re-order columns so UID is first
db$UID <- 1:nrow(db)
db <- db[,c(ncol(db),1:(ncol(db)-1))]
rownames(db) <- db$UID

#--------------

#Remove samples and taxa with no observations
#-Identify and remove samples with no diatom counts (i.e., rowSum = 0)
db[rowSums(db[,9:ncol(db)]) == 0,][,1:9]
db <- db[rowSums(db[,9:ncol(db)]) > 0,]

#-Identify and remove taxa with no counts across all samples (i.e., colSum = 0)
colnames(db[,9:ncol(db)][,colSums(db[,9:ncol(db)]) == 0])
db <- db[, c(1:8, (which(colSums(db[,9:ncol(db)]) > 0)+8))]

#--------------

#Convert counts to relative abundance, dividing each count-per-row by the rowSum and multiplying by 100
db_rel <- db
db_rel[,9:ncol(db)] <- db_rel[,9:ncol(db)] / rowSums(db_rel[,9:ncol(db)]) * 100

#--------------

#Extract relative abundance table in long format for input into MMI Shiny application
#-Convert db from wide to long format, retaining only UID, taxon name, and relative abundance for each taxon-by-UID
db_long <- reshape2::melt(db_rel[,c(1,9:ncol(db_rel))], id.vars = c("UID"))

#-Rename columns to match MMI Shiny input
colnames(db_long) <- c("SAMPLEID", "TAXON", "COUNT")

#-Select only the records for which a taxon was recorded (i.e., relative abundance > 0)
db_long <- db_long[which(db_long$COUNT > 0),]

#--------------

#Download relative abundance of counts in long format (db_long) for input into MMI Shiny application (https://diatom-mmi-demo.shinyapps.io/Diatom_MMI_V24/)
#-Save 'db_long' as "MMI_Input_Step_2_Sample_Taxonomy_Counts.csv"
#--Upload this file into the MMI Shiny Application during 'Step 2: Upload Taxonomic Data'
#--As of 05/2021, the taxonomy list is as updated as is possible. 
#--If no updates in CDPHE or MMI taxonomy since 05/2021, proceed with 'Step 3 MMI Calculations' without further harmonizing data (following Step 1 - processing outlined below)
#write.csv(db_long, "MMI_Input_Step_2_Sample_Taxonomy_Counts.csv", row.names = FALSE)

#--------------

print(noquote("Processing of diatom data prior to input into MMI is complete!"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#####Process sample information prior to input into MMI

#--------------

#Extract site metadata, add COMID and drainage area, and re-format to match input requirements for MMI Shiny application

#-Extract metadata from WQCD Algal Database
meta <- db[,1:8]

#-Load station information with COMIDs (National Hydrography Dataset unique stream segment identifier)
#--COMIDs were assigned manually and verified as part of this project
comid <- read.csv("CDPHE_COMID.csv", header = TRUE, strip.white = TRUE)

#-Load station information with drainage area
#--Drainage area was calculated previously by CDPHE, with some areas calculated as part of this project
drainage <- read.csv("CDPHE_Drainage.csv", header = TRUE, strip.white = TRUE)

#-Add COMID and drainage area to metadata based on StationID
meta$COMID <- comid[match(meta$StationID, comid$StationID), "COMID"]
meta$catchment_km2 <- drainage[match(meta$StationID, drainage$StationID), "Catchment_Area_SQ_KM"]

#-Retain only COMID, UID, and drainage area, re-order, and re-name for input into MMI Shiny application
site <- meta[,c(9,1,10)]
colnames(site) <- c("COMID", "SAMPLEID", "AREA")

rm(comid, drainage)

#Download site information for input into MMI Shiny application (https://diatom-mmi-demo.shinyapps.io/Diatom_MMI_V24/)
#-Save 'site' as "MMI_Input_Step_1_Site_Information.csv"
#--Upload this file into the MMI Shiny Application during 'Step 1: Upload Site Information'
#--As of 05/2021, the COMIDs and drainage areas are updated 
#--If no updates in StationIDs and/or locations since 05/2021, proceed with 'Step 3 MMI Calculations' without further harmonizing data
#write.csv(site, "MMI_Input_Step_1_Site_Information.csv", row.names = FALSE)

#--------------

print(noquote("Processing of site information prior to input into MMI is complete!"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#####Computing the MMI

#--------------

#All MMI computations are performed in the MMI Shiny application. Refer to report text for detailed instructions
#To proceed with script below for MMI interpretation and environmental analysis:
#-Download "PLAINS & LOWLANDS MMI OUTPUT" and "WEST MMI OUTPUT"
#-Move MMI output files to working directory

print(noquote("DO NOT PROCEED UNLESS MMI OUTPUT FILES HAVE BEEN DOWNLOADED FROM SHINY AND MOVED TO WORKING DIRECTORY"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#####Process environmental data prior to analysis of MMI scores and diatom community patterns

#--------------

#Align water chemistry data with count data

#-Extract UID and relative abundance data from db_rel
count <- db_rel[,c(1,9:ncol(db_rel))]

#-Create unique sample identifiers by merging Station ID and collection date to be able to match to chemistry dataset that lacks UIDs
meta$Names <- paste(gsub("/", ".", meta$`Sample Date`), meta$StationID, sep = "_")

#-Remove duplicate samples
meta <- meta[-grep("Duplicate", meta$`Sample Type`), ] #First remove all samples with SampleType designated as "Duplicate"
meta <- meta[!duplicated(meta$Names), ] #Next, remove duplicate samples with no "Duplicate" designated

#-Load chemistry data
print(noquote("Note that in WQCD Stressor ID Pairings any value with red font was turned to blank in Excel prior to saving as 'WQCD Stressor ID Pairings - 070720_rm.csv' and uploading"))
chem <- read.csv("WQCD Stressor ID Pairings - 070720_rm.csv", header = TRUE, strip.white = TRUE, na.strings = c("", NA))

#-Create unique sample identifiers by merging Station ID and collection date to match to metadata
chem$Names <- paste(gsub("/", ".", chem$CollDate), chem$StationID, sep = "_")

#-Remove duplicate 'chem' samples
chem <- chem[!duplicated(chem$Names), ]

#-Match UID to chemistry data
chem$UID <- meta[match(chem$Names, meta$Names), "UID"]

#-Move UID and Names columns to before chemistry data
chem <- chem[, c(ncol(chem), (ncol(chem)-1), 1:(ncol(chem)-2))]

#-Retain full chemistry dataset with UIDs and duplicates removed (used for map making)
chem_full <- chem

#-Retain only chemistry samples with count data, and only count samples with chemistry data
chem <-  chem[(chem$UID %in% count$UID), ]
count <-  count[(count$UID %in% chem$UID), ]

#-Re-name 'chem' rownames to UID
rownames(chem) <- chem$UID

print(noquote(paste(nrow(chem), " of ", nrow(db), " diatom count samples contain environmental data", sep = "")))

#--------------

#Process chemistry data for BDL and missing values
#-Replace all instances of "broke"
chem[chem == "broke"] <- NA

#-Count percent of BDLs per variable
bdl_vec <- vector()
for(a in colnames(chem[,9:ncol(chem)])){
  bdl_vec[[a]] <- length(grep("<", chem[[a]])) / nrow(chem) * 100
}
bdl_vec

#-Replace "<x" with =x/2 (half of detection limit as specified by "<x")
for(j in 9:ncol(chem)){
  for(i in 1:nrow(chem)){
    chem[i,j] <- ifelse(grepl("<", chem[i,j]),
                        as.numeric(sub('.', '', chem[i,j])) / 2,
                        as.numeric(chem[i,j]))
  }
}

#-Change all variable columns to numeric
chem[,9:ncol(chem)] <- data.frame(lapply(chem[,9:ncol(chem)], function(x) as.numeric(as.character(x))))

#-Count percent of NAs per column
na_vec <- vector()
for(a in colnames(chem[,9:ncol(chem)])){
  na_vec[[a]] <- sum(is.na(chem[[a]])) / nrow(chem) * 100
}
na_vec

#Create dataframe of %BDL and %NA per variable of the 351 samples with count and chem data
missing_chem <- as.data.frame(cbind(bdl_vec, na_vec))
rownames(missing_chem) <- gsub("\\.L", "/L", rownames(missing_chem))
rownames(missing_chem) <- gsub("\\.", " ", rownames(missing_chem))
rownames(missing_chem) <- gsub(" cm", "/cm", rownames(missing_chem))
rownames(missing_chem) <- gsub("  ", " ", rownames(missing_chem))
colnames(missing_chem)[1] <- "BDL (%)"
colnames(missing_chem)[2] <- "NA (%)"
missing_chem

#--------------

#Transform each environmental variable in a way that minimizes skewness and kurtosis
#-Create copy of 'chem', retaining only environmental variable columns
chem_trans <- chem[, 9:ncol(chem)]

#-Create empty dataframe for recording which transformations were performed
chem_trans_meta <- data.frame(transformation = rep(NA, ncol(chem_trans)),
                              skewness = rep(NA, ncol(chem_trans)),
                              kurtosis = rep(NA, ncol(chem_trans)))
rownames(chem_trans_meta) <- colnames(chem_trans)

#-Perform transformations from list of possible options and populate list of transformations that were performed
library(moments)
for(a in colnames(chem_trans)){
  df <- data.frame(transformation = c("quart", "cube", "square", "no transformation", "square root", "cube root", "fourth root", "log(x)", "reciprocal root", "reciprocal", "reciprocal square"),
                   formula = c("chem_trans[[a]]^4", "chem_trans[[a]]^3", "chem_trans[[a]]^2", "chem_trans[[a]]", "chem_trans[[a]]^0.5", "chem_trans[[a]]^(1/3)", "chem_trans[[a]]^0.25", "log(chem_trans[[a]])", "-1 / chem_trans[[a]]^0.5", "-1 / chem_trans[[a]]", "-1 / chem_trans[[a]]^2"),
                   skewness = rep(NA, 11),
                   kurtosis = rep(NA, 11),
                   skew_kurt_sum = rep(NA, 11))
  for(i in 1:nrow(df)){
    df[i,3] <- skewness(eval(parse(text = paste0(df[i,2]))), na.rm = TRUE)
    df[i,4] <- kurtosis(eval(parse(text = paste0(df[i,2]))), na.rm = TRUE)
    df[i,5] <- abs(skewness(eval(parse(text = paste0(df[i,2]))), na.rm = TRUE)) + abs(kurtosis(eval(parse(text = paste0(df[i,2]))), na.rm = TRUE))
  }
  df <- df[order(df$skew_kurt_sum),]
  df <- as.matrix(head(df,1))
  chem_trans[[a]] <- eval(parse(text = paste0(df[1,2])))
  chem_trans_meta[a,1] <- df[1,1]
  chem_trans_meta[a,2] <- df[1,3]
  chem_trans_meta[a,3] <- df[1,4]
}

#-Examine histograms for each untransformed and transformed variable to manually check the reduction of skewness and kurtosis
for(a in colnames(chem_trans)){
  hist(chem[[a]], main = paste(a, ": untransformed", sep = ""))
  hist(chem_trans[[a]], main = paste(a, ": transformed", sep = ""))
}

#-Save transformations performed on each variable
#write.csv(chem_trans_meta, "Variable_Transformations.csv")

#--------------

#Select variable to retain for analysis based on correlations, missing values, ecological relevance, and inclusivity
#-Calculate Spearman correlations for complete cases of variables
corrs <- cor(chem_trans[,1:ncol(chem_trans)], method = c("spearman"), use = "complete.obs")
rownames(corrs) <- gsub("\\.cm.", "/cm)", rownames(corrs))
rownames(corrs) <- gsub("\\.L.", "/L)", rownames(corrs))
rownames(corrs) <- gsub("\\.\\.", " (", rownames(corrs))
rownames(corrs) <- gsub("\\.", "_", rownames(corrs))
colnames(corrs) <- gsub("\\.cm.", "/cm)", colnames(corrs))
colnames(corrs) <- gsub("\\.L.", "/L)", colnames(corrs))
colnames(corrs) <- gsub("\\.\\.", " (", colnames(corrs))
colnames(corrs) <- gsub("\\.", "_", colnames(corrs))

#-Optional: Remove correlations < |0.5| for easier visualization
#corrs[corrs < abs(0.5)] <- 0

#-Visualize correlations in correlation plot
pdf(file = "CDPHE_Water_Quality_Correlations.pdf")
corrplot::corrplot(corrs, tl.cex = 0.5, type = "upper", tl.col = "black", 
                   col = colorRampPalette(c("gray50","white","gray49"))(10),
                   addCoef.col = "black", number.cex = 0.5, diag = FALSE)
dev.off()

#-Based on correlations (|x| >= 0.7) and missing values, select the following variables for further analysis:
#--PH
#--COND (correlated with alkalinity, chloride, selenium)
#--COPPER
#--IRON_DIS
#--AMMONIA
#--PHOSPHORUS (correlated with NO5, temperature, chloride)
#--ZINC
#----NO5 removed b/c 0.71 correlation with PHOSPHORUS
chem_trans <- chem_trans[, c(1:2,7,10,13,17,20)]

#Remove samples with NA
chem_trans <- chem_trans[complete.cases(chem_trans[, 1:ncol(chem_trans)]), ]

#--------------

#Retain samples with both count and filtered chemistry data
#-Retain only count samples with environmental data, and environmental samples with count data
count <-  count[(count$UID %in% rownames(chem_trans)), ]
chem_trans <-  chem_trans[(rownames(chem_trans) %in% count$UID), ]

#Clean up and re-order 'chem_trans' and 'count' dataframes to match
chem_trans <- chem_trans[gtools::mixedorder(rownames(chem_trans)), ]
count <- count[gtools::mixedorder(rownames(count)), ]
count$UID <- NULL

#Create separate df for metadata and remove metadata from chem_trans
meta <- meta[rownames(meta) %in% rownames(chem_trans), ]
meta <- meta[gtools::mixedorder(rownames(meta)), ]

#Remove objects no longer needed
rm(corrs, df, a, bdl_vec, i, j, na_vec)

print(noquote(paste(nrow(chem_trans), " of ", nrow(db), " diatom count samples with no missing data for selected environmental variables", sep = "")))

print(noquote("Processing of environmental data prior to data analysis is complete!"))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#####Analyze MMI scores and diatom community patterns using environmental data

#--------------

#Construct non-metric multidimensional scaling (NMDS) ordination for relative abundance of diatom counts
#-Remove outlier samples (determined by first running an NMDS with all samples and identifying outlier points)
rows_rm <- c("295", "343", "344")
count_rm <- count[!(rownames(count) %in% rows_rm), ]
chem_trans_rm <- chem_trans[!(rownames(chem_trans) %in% rows_rm), ]
meta_rm <- meta[!(rownames(meta) %in% rows_rm), ]

#-Run NMDS
set.seed(10) #Ensure the same NMDS is constructed every time
nmds_count <- vegan::metaMDS(count_rm, distance = "bray", trymax = 99)

#-Extract site scores
site_scores_count <- as.data.frame(vegan::scores(nmds_count, "sites"))

#-Extract species scores
species_scores_count <- as.data.frame(vegan::scores(nmds_count, "species"))
species_scores_count <- species_scores_count[complete.cases(species_scores_count), ]
species_scores_count$taxon <- rownames(species_scores_count)

#--Correlate species relative abundance with NMDS points
corrs <- as.data.frame(cor(count_rm, nmds_count$points, use = "complete.obs", method = "spearman")) #Ignore any warnings of 'the standard deviation is zero'
corrs <- corrs[complete.cases(corrs), ]

#--Retain species scores if species relative abundance is correlated at |??| >= 0.3 with one of the two NMDS axes (for visualization only)
species_retain <- corrs[abs(corrs$MDS1) >= 0.3 | abs(corrs$MDS2) >= 0.3, ]
species_scores_count <- species_scores_count[(species_scores_count$taxon %in% rownames(species_retain)), ]

#-Generate environmental vectors
set.seed(10)
env_count <- vegan::envfit(nmds_count, chem_trans_rm, perm = 99)
env_count$vectors #To see R2 and significance of vectors
env_seg_count <- as.data.frame(env_count$vectors$arrows*sqrt(env_count$vectors$r))
env_seg_count$chem <- c("pH", "Conductivity", "Cu", "Fe", "NH3", "TP", "Zn")
env_seg_count$r2 <- env_count$vectors$r
env_seg_count$pvals <- env_count$vectors$pvals

#--Retain only environmental vectors that are significant or explanatory (R2 >= 0.15) (for visualization only)
env_seg_count <- env_seg_count[env_seg_count$pvals < 0.05, ]
env_seg_count <- env_seg_count[env_seg_count$r2 >= 0.15, ]
env_seg_count$chem <- paste(env_seg_count$chem, " (", round(env_seg_count$r2, digits = 2), ")", sep = "")

#--Calculate points so env labels are at end of arrows
library(tidyverse)
library(patchwork)
#---Radial shift function
rshift = function(r, theta, a=0.03, b=0.07) {
  r + a + b*abs(cos(theta))
}
#---Calculate shift
env_seg_count <- env_seg_count %>% 
  mutate(r = sqrt(NMDS1^2 + NMDS2^2),
         theta = atan2(NMDS2,NMDS1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

#-Plot NMDS ordination
library(ggplot2)
library(ggrepel)
plot_count <- ggplot() + 
  geom_point(data = site_scores_count, aes(x = NMDS1, y = NMDS2), shape = 21, size = 2, color = "gray30", fill = "gray70", alpha = 0.15) + # add the point markers
  geom_segment(data = env_seg_count, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2), #Draw environmental vector arrows
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"), size = 1, color = "gray30") + 
  geom_text_repel(data = species_scores_count, aes(x = NMDS1, y = NMDS2, label = taxon), size = 2.5, color = "gray30", fontface = "italic") + #Write species_scores_counts of species correlated with axes
  geom_text(data = env_seg_count, aes(x = xnew, y = ynew, label = chem), size = 3, color = "gray30", fontface = "bold") + #Write environmental variable labels
  annotate(geom = 'text', label = paste("Stress = ", round(nmds_count$stress, digits = 2)), #Add stress value to top right
           x = Inf, y = Inf, hjust = "right", vjust = "top") +
  scale_x_continuous(expand = c(.1, .1)) +
  scale_y_continuous(expand = c(.1, .1)) +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_blank(), # remove x-axis labels
        axis.title.y = element_blank(), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        aspect.ratio = 1/1)
plot(plot_count)
ggsave("nmds_count.png",
       plot = last_plot(),
       width = 3.25,
       height = 3.25,
       units = c("in"),
       dpi = 1200)

rm(chem_trans_rm, corrs, count_rm, env_count, env_seg_count, meta_rm, nmds_count, rows_rm, site_scores_count, species_retain, species_scores_count, rshift)

#--------------

#Construct NMDS ordination for scores of metrics in MMIs
#-Plains and Lowlands (ELO)
#--Load metric data from ELO
#---This file is the output from "Step 3: MMI Calculations" > Download "Plains & Lowlands MMI Output" in the MMI Shiny application
#----Move this .csv file from your Downloads folder to your working directory prior to continuing
elo <- read.csv("ELO_MMI.csv", header = TRUE, strip.white = TRUE)

#--Edit rownames
rownames(elo) <- elo$SAMPLEID

#--Invert metrics that were re-scaled during MMI computation
#---NOTE: These metrics may need to change if component metrics change in future versions of the MMI
#----As a result, any command that refers to the MMI output files should be checked for consistency 
#-----E.g., removing non-metric columns may be columns 1 and 7 now, but could be 1 and 8 if the MMI is updated with 6 component metrics)
elo$SALINITY_3.pt.res <- 10 - elo$SALINITY_3.pt.res
elo$CENTRIC.pt <- 10 - elo$CENTRIC.pt

#--Retain full dataframe with inverted metrics
elo_full <- elo

#--Remove non-metric columns
elo <- elo[, -c(1,7)]

#--Retain only ELO samples with environmental data and re-order
elo <-  elo[(rownames(elo) %in% rownames(chem_trans)), ]
elo <- elo[gtools::mixedorder(rownames(elo)), ]

#--Run NMDS
set.seed(10)
nmds_elo <- vegan::metaMDS(elo, distance = "bray", trymax = 99)

#--Extract site scores
site_scores_elo <- as.data.frame(vegan::scores(nmds_elo, "sites"))

#--Extract species scores
species_scores_elo <- as.data.frame(vegan::scores(nmds_elo, "species"))
species_scores_elo <- species_scores_elo[complete.cases(species_scores_elo), ]
species_scores_elo$taxon <- rownames(species_scores_elo)

#--Generate environmental vectors
elo_chem <- chem_trans[(rownames(chem_trans) %in% rownames(elo)), ]
elo_chem <- elo_chem[gtools::mixedorder(rownames(elo_chem)), ]
set.seed(10)
env_elo <- vegan::envfit(nmds_elo, elo_chem, perm = 100)
head(env_elo)
env_seg_elo <- as.data.frame(env_elo$vectors$arrows*sqrt(env_elo$vectors$r))
env_seg_elo$chem <- c("pH", "Conductivity", "Cu", "Fe", "NH3", "TP", "Zn")
env_seg_elo$r2 <- env_elo$vectors$r
env_seg_elo$pvals <- env_elo$vectors$pvals

#---Retain only environmental vectors that are significant or explanatory (R2 >= 0.15) (for visualization only)
env_seg_elo <- env_seg_elo[env_seg_elo$pvals < 0.05, ]
env_seg_elo <- env_seg_elo[env_seg_elo$r2 >= 0.15, ]
env_seg_elo$chem <- paste(env_seg_elo$chem, " (", round(env_seg_elo$r2, digits = 2), ")", sep = "")

#---Calculate points so env labels are at end of arrows (Ref: https://stackoverflow.com/questions/64935396/how-do-i-shift-the-geom-text-labels-to-after-a-geom-segment-arrow-in-ggplot2)
library(tidyverse)
library(patchwork)
#----Radial shift function
rshift = function(r, theta, a=0.03, b=0.07) {
  r + a + b*abs(cos(theta))
}
#----Calculate shift
env_seg_elo <- env_seg_elo %>% 
  mutate(r = sqrt(NMDS1^2 + NMDS2^2),
         theta = atan2(NMDS2,NMDS1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

#--Plot NMDS ordination
library(ggplot2)
library(ggrepel)
plot_elo <- ggplot() + 
  geom_point(data = site_scores_elo, aes(x = NMDS1, y = NMDS2), shape = 21, size = 2, color = "gray30", fill = "gray70", alpha = 0.15) +
  geom_segment(data = env_seg_elo, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.15, "cm")), size = 1, color = "gray30") + 
  geom_text_repel(data = species_scores_elo, aes(x = NMDS1, y = NMDS2, label = taxon), size = 2.5, color = "gray30") +
  geom_text(data = env_seg_elo, aes(x = xnew, y = ynew, label = chem), size = 3, color = "gray30", fontface = "bold") +
  annotate(geom = 'text', label = paste("Stress = ", round(nmds_elo$stress, digits = 2)), 
           x = Inf, y = Inf, hjust = "right", vjust = "top") +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_blank(), # remove x-axis labels
        axis.title.y = element_blank(), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        aspect.ratio = 1/1)
plot(plot_elo)
ggsave("nmds_elo.png",
       plot = last_plot(),
       width = 3.25,
       height = 3.25,
       units = c("in"),
       dpi = 1200)

rm(env_elo, env_seg_elo, nmds_elo, site_scores_elo, species_scores_elo, rshift)

#---

#-West (WST)
#--Load metric data from WST
#---This file is the output from "Step 3: MMI Calculations" > Download "West MMI Output" in the MMI Shiny application
#----Move this .csv file from your Downloads folder to your working directory prior to continuing
wst <- read.csv("WST_MMI.csv", header = TRUE, strip.white = TRUE)

#--Edit rownames
rownames(wst) <- wst$SAMPLEID

#--Invert metrics that were re-scaled during MMI computation
wst$BC_4.r.res <- 10 - wst$BC_4.r.res
wst$HIGH_N.pa.res <- 10 - wst$HIGH_N.pa.res
wst$SALINITY_3.pt <- 10 - wst$SALINITY_3.pt

#--Retain full dataframe with inverted metrics
wst_full <- wst

#--Remove non-metric columns
wst <- wst[, -c(1,8)]

#--Retain only wst samples with environmental data and re-order
wst <-  wst[(rownames(wst) %in% rownames(chem_trans)), ]
wst <- wst[gtools::mixedorder(rownames(wst)), ]

#--Run NMDS
set.seed(10)
nmds_wst <- vegan::metaMDS(wst, distance = "euclidean", trymax = 99)

#--Extract site scores
site_scores_wst <- as.data.frame(vegan::scores(nmds_wst, "sites"))

#--Extract species scores
species_scores_wst <- as.data.frame(vegan::scores(nmds_wst, "species"))
species_scores_wst <- species_scores_wst[complete.cases(species_scores_wst), ]
species_scores_wst$taxon <- rownames(species_scores_wst)

#--Generate environmental vectors
wst_chem <- chem_trans[(rownames(chem_trans) %in% rownames(wst)), ]
wst_chem <- wst_chem[gtools::mixedorder(rownames(wst_chem)), ]
set.seed(10)
env_wst <- vegan::envfit(nmds_wst, wst_chem, perm = 100)
head(env_wst)
env_seg_wst <- as.data.frame(env_wst$vectors$arrows*sqrt(env_wst$vectors$r))
env_seg_wst$chem <- c("pH", "Conductivity", "Cu", "Fe", "NH3", "TP", "Zn")
env_seg_wst$r2 <- env_wst$vectors$r
env_seg_wst$pvals <- env_wst$vectors$pvals

#---Remove environmental vectors if not significant (for visualization only)
env_seg_wst <- env_seg_wst[env_seg_wst$pvals < 0.05, ]
env_seg_wst <- env_seg_wst[env_seg_wst$r2 >= 0.15, ]
env_seg_wst$chem <- paste(env_seg_wst$chem, " (", round(env_seg_wst$r2, digits = 2), ")", sep = "")

#---Calculate points so env labels are at end of arrows
library(tidyverse)
library(patchwork)
#----Radial shift function
rshift = function(r, theta, a=0.03, b=0.07) {
  r + a + b*abs(cos(theta))
}
#----Calculate shift
env_seg_wst <- env_seg_wst %>% 
  mutate(r = sqrt(NMDS1^2 + NMDS2^2),
         theta = atan2(NMDS2,NMDS1),
         rnew = rshift(r, theta),
         xnew = rnew*cos(theta),
         ynew = rnew*sin(theta))

#--Plot NMDS ordination
library(ggplot2)
library(ggrepel)
plot_wst <- ggplot() + 
  geom_point(data = site_scores_wst, aes(x = NMDS1, y = NMDS2), shape = 21, size = 2, color = "gray30", fill = "gray70", alpha = 0.15) +
  geom_segment(data = env_seg_wst, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0.15, "cm")), size = 1, color = "gray30") + 
  geom_text_repel(data = species_scores_wst, aes(x = NMDS1, y = NMDS2, label = taxon), size = 2.5, color = "gray30") +
  geom_text(data = env_seg_wst, aes(x = xnew, y = ynew, label = chem), size = 3, color = "gray30", fontface = "bold") +
  annotate(geom = 'text', label = paste("Stress = ", round(nmds_wst$stress, digits = 2)), 
           x = Inf, y = Inf, hjust = "right", vjust = "top") +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_blank(), # remove x-axis labels
        axis.title.y = element_blank(), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        plot.margin=unit(c(0,0,0,0),"mm"),
        aspect.ratio = 1/1)
plot(plot_wst)
ggsave("nmds_wst.png",
       plot = last_plot(),
       width = 3.25,
       height = 3.25,
       units = c("in"),
       dpi = 1200)

rm(env_wst, env_seg_wst, nmds_wst, site_scores_wst, species_scores_wst, rshift)

egg::ggarrange(plot_count, plot_wst, plot_elo, ncol = 3)
ggsave("test.png",
       plot = last_plot(),
       width = 9.75,
       height = 3.25,
       units = c("in"),
       dpi = 1200)

#--------------

#Perform linear regressions and construct generalized additive models (GAMs) for MMI scores
#-Plains and Lowlands (ELO)
#--Create dataframe of MMI scores and environmental data
elo_gam <- elo_full
elo_gam <-  elo_gam[(rownames(elo_gam) %in% rownames(chem_trans)), ]
elo_gam <- elo_gam[gtools::mixedorder(rownames(elo_gam)), ]
elo_global <- cbind(elo_gam, elo_chem)

#--Perform simple linear regressions for each environmental variable
elo_lm <- data.frame(Response = NA, Explanatory = NA, R2 = NA, p = NA)
for(i in 1:ncol(elo_chem)){
  a <- colnames(elo_chem)[i]
  temp_lm <- lm(elo_global$MMI ~ elo_global[[a]])
  elo_lm[i,1] <- "MMI"
  elo_lm[i,2] <- a
  elo_lm[i,3] <- summary(temp_lm)$adj.r.squared
  elo_lm[i,4] <- summary(temp_lm)$coefficients[2,4]
  r2 <- round(summary(temp_lm)$adj.r.squared, digits = 2)
  p <- round(summary(temp_lm)$coefficients[2,4], digits = 2)
  print(ggResidpanel::resid_panel(temp_lm))
  print(ggplot() +
          geom_point(aes(elo_global[[a]], elo_global$MMI), shape = 21, size = 2, color = "gray30", fill = "gray70", alpha = 0.50) +
          geom_smooth(aes(elo_global[[a]], elo_global$MMI), method = lm, se = TRUE, color = "gray30", fill = "gray", alpha = 0.50) +
          #labs(x = paste(a, "(transformed)"), y = "MMI score") +
          annotate(geom = "label", label = deparse(bquote(atop(paste(R^2 ~ "=" ~ .(r2)), italic(.("p") ~ "=" ~ .(p))))), parse = TRUE,
                   label.r = unit(0, "lines"), size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1, label.size = 0.25, label.padding = unit(0.075, "lines")) +
          scale_x_continuous(expand = c(0, 0)) +
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank()))
  ggsave(paste("elo_lm_", a, ".png", sep = ""),
         plot = last_plot(),
         width = 2.4,
         height = 2.4,
         units = c("in"),
         dpi = 1200)
}
elo_lm

#--Perform multiple linear regression for all environmental variables
elo_mlm <- lm(MMI ~ ., data = elo_global[,7:14])
print(ggResidpanel::resid_panel(elo_mlm))
summary(elo_mlm)

#--Construct all possible GAMs with 1-2 explanatory variables
library(mgcv)

#---Create all 1-2 variable combinations using all variables
elo_combos <- list()
elo_combvec <- vector()

for(a in colnames(elo_chem)){
  value <- paste0("s(", a, ", k = 5)")
  elo_combvec[[a]] <- value
}
elo_combos$one <- as.matrix(t(combn(elo_combvec, 1)))
for(i in 1:nrow(elo_combos$one)){
  value <- paste(elo_combos$one[i,], collapse = " + ")
  elo_combos$one[i,1] <- value
}
elo_combos$one <- elo_combos$one[,-2]
elo_combos$two <- as.matrix(t(combn(elo_combvec, 2)))
for(i in 1:nrow(elo_combos$two)){
  value <- paste(elo_combos$two[i,], collapse = " + ")
  elo_combos$two[i,1] <- value
}
elo_combos$two <- elo_combos$two[,-2]
elo_combos_merged <- as.matrix(do.call(c, elo_combos))

rm(elo_combvec, value)

#---Add response variable name to variable combinations
elo_combos_merged <- as.matrix(paste("MMI ~", elo_combos_merged))

#---Run GAM all subsets
n <- nrow(elo_combos_merged)
elo_mods <- data.frame(formula = rep(NA, n),
                       max.pval = rep(NA, n),
                       dev.exp = rep(NA, n),
                       concurvity.worst = rep(NA, n),
                       concurvity.obs = rep(NA, n),
                       concurvity.est = rep(NA, n),
                       AIC = rep(NA, n))
for(i in 1:nrow(elo_combos_merged)){
  mod <- gam(as.formula(paste(elo_combos_merged[i,1])),
             data = elo_global,
             fx = FALSE, #Set terms to penalized regression splines
             select = TRUE, #Give each smooth an extra penalty on its fixed effect component (penalize null space of each smooth)
             gamma = 1.4, #Increase penalty on each model degree of freedom to suppress overfitting
             family = tw, #set distribution family (Tweedie)
             method = "REML", #set smoothing parameter estimation
             na.action = "na.fail")
  con <- concurvity(mod)
  elo_mods[i, 1] <- paste0(summary(mod)$formula[2],summary(mod)$formula[1],summary(mod)$formula[3])
  elo_mods[i, 2] <- max(summary(mod)$s.pv)
  elo_mods[i, 3] <- summary(mod)$dev.expl
  elo_mods[i, 4] <- max(con[1, -1])
  elo_mods[i, 5] <- max(con[2, -1])
  elo_mods[i, 6] <- max(con[3, -1])
  elo_mods[i, 7] <- AIC(mod)
}

elo_stored <- elo_mods #Retain all subsets models for reference

#---Filter out nonsignificant, low deviance explained, and low performing GAMs
elo_mods <- elo_mods[elo_mods$max.pval <= 0.05, ]
elo_mods <- elo_mods[elo_mods$dev.exp >= 0.15, ]
elo_mods <- elo_mods[elo_mods$concurvity.est < 0.70, ]
elo_mods$concurvity.worst <- NULL
elo_mods$concurvity.obs <- NULL

#---Determine models with significant smooth terms or non-convergence in gam.check (also look at residual plots)
elo_mods$min_pval_check <- NA
elo_mods$convergence <- NA
for(i in 1:nrow(elo_mods)){
  mod <- gam(as.formula(paste(elo_mods[i,1])),
             data = elo_global,
             fx = FALSE, 
             select = TRUE, 
             gamma = 1.4, 
             family = tw, 
             method = "REML", 
             na.action = "na.fail")
  set.seed(10)
  elo_mods[i,6] <- min(k.check(mod)[,4]) #Minimum p-value of smooth terms in gam.check
  set.seed(10)
  elo_mods[i,7] <- capture.output(gam.check(mod))[3] #Whether model converged
}

#---Calculate delta AIC and Akaike weights and likelihoods
library(qpcR)
elo_mods <- elo_mods[elo_mods$min_pval_check > 0.05, ] #Remove models with significant smooth terms in gam.check (i.e., oversmoothed)
elo_mods <- elo_mods[grep("full convergence", elo_mods$convergence), ] #Remove re-run models that did not converge
elo_mods <- elo_mods[elo_mods$max.pval < 0.05, ]
elo_mods <- elo_mods[elo_mods$dev.exp >= 0.15, ]
elo_mods$delta.AIC <- akaike.weights(elo_mods$AIC)$deltaAIC
elo_mods$akaike.likelihoods <- akaike.weights(elo_mods$AIC)$rel.LL
elo_mods$akaike.weights <- akaike.weights(elo_mods$AIC)$weights
elo_mods <- elo_mods[elo_mods$delta.AIC < 2, ] #Retain only those models with delta AIC < 2
elo_mods <- elo_mods[order(-elo_mods$akaike.weights), ] #Order models by Akaike weight

#---Check GAM residuals to validate model
for(i in 1:nrow(elo_mods)){
  mod <- gam(as.formula(paste(elo_mods[i,1])),
             data = elo_global,
             fx = FALSE,
             select = TRUE,
             gamma = 1.4,
             family = tw,
             method = "REML",
             na.action = "na.fail")
  gam.check(mod)
}

#---Un-transform explanatory variable values for partial effects plots
df <- data.frame(transformation = c("quart", "cube", "square", "no transformation", "square root", "cube root", "fourth root", "log(x + 1)", "reciprocal root + 1", "reciprocal + 1", "reciprocal square +1", "log(x)", "reciprocal root", "reciprocal", "reciprocal square"),
                 formula = c("elo_viz[[i]][1]^0.25", "elo_viz[[i]][1]^(1/3)", "elo_viz[[i]][1]^0.5", "elo_viz[[i]][1]", "elo_viz[[i]][1]^2", "elo_viz[[i]][1]^3", "elo_viz[[i]][1]^4", "exp(elo_viz[[i]][1]) - 1", "(-1 / elo_viz[[i]][1])^2 - 1", "(-1 / elo_viz[[i]][1]) - 1", "(-1 / elo_viz[[i]][1])^0.5 - 1","exp(elo_viz[[i]][1])", "(-1 / elo_viz[[i]][1])^2", "-1 / elo_viz[[i]][1]", "(-1 / elo_viz[[i]][1])^0.5"),
                 blank = rep(NA, 15))
chem_trans_meta$elo_inverse <- lapply(chem_trans_meta$transformation, function(x) df$formula[match(x, df$transformation)])

#---Construct partial effects plots for filtered GAMs
#----Re-run model
mod <- gam(as.formula(paste(elo_mods[1,1])),
           data = elo_global,
           fx = FALSE,
           select = TRUE,
           gamma = 1.4,
           family = tw,
           method = "REML",
           na.action = "na.fail")

#----Extract fitted values and residuals
library(mgcViz)
elo_viz <- list()
elo_viz_res <- list()
for(i in 1:length(mod$var.summary)){
  var <- names(mod$var.summary)[[i]] 
  o <- plot(sm(getViz(mod), i))
  elo_viz[[var]] <- o$data$fit
  elo_viz_res[[var]] <- o$data$res
}

#----Process and plot fitted values and residuals
for(i in names(elo_viz)){
  a <- i
  elo_viz[[i]][[a]] <- elo_viz[[i]]$x
  elo_viz[[i]][[a]] <- eval(parse(text = paste0(as.data.frame(t(chem_trans_meta))[4,a])))
  elo_viz[[i]] <- as.data.frame(as.matrix(elo_viz[[i]]))
  
  elo_viz_res[[i]][[a]] <- elo_viz_res[[i]]$x
  elo_viz_res[[i]][[a]] <- eval(parse(text = gsub("viz", "viz_res", paste0(as.data.frame(t(chem_trans_meta))[4,a]))))
  elo_viz_res[[i]] <- as.data.frame(as.matrix(elo_viz_res[[i]]))
  
  plot(ggplot(elo_viz[[i]], aes(x = elo_viz[[i]][[a]], y = y)) +
         geom_ribbon(aes(ymin = y - se, ymax = y + se), fill = "gray", color = NA, alpha = 0.15) + #Plot confidence intervals
         #geom_point(data = elo_viz_res[[i]], aes(x = elo_viz_res[[i]][[a]], y = y), shape = 21, size = 2, color = "gray50", fill = "gray70") + #Optional: plot residuals that fall within fitted line
         geom_rug(data = elo_viz_res[[i]], aes(x = elo_viz_res[[i]][[i]], y = NULL), color = "gray30") + #Plot rugs along x-axis to represent distribution of all residuals
         geom_line(size = 1, color = "gray30") + #Plot fitted line
         #labs(title = paste("Model deviance explained = ", round(summary(mod)$dev.expl, digits = 2), sep = ""), x = a, y = "Partial effect") +
         scale_x_continuous(limits = c(min(elo_viz[[i]][[a]]), max(elo_viz[[i]][[a]])), expand = c(0, 0)) +
         scale_y_continuous(limits = c(min(elo_viz[[i]]$y - elo_viz[[i]]$se), max(elo_viz[[i]]$y + elo_viz[[i]]$se)), expand = c(0, 0)) + #Only display residuals within the y-axis bounds of the fit
         theme_bw()+
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank()))
  
  ggsave(paste("elo_gam_", i, ".png", sep = ""),
         plot = last_plot(),
         width = 2.4,
         height = 2.4,
         units = c("in"),
         dpi = 1200)
}

print(noquote("If 'elo_mods' has 0 observations, no GAMs are significant and explanatory. Disregard any plots or relationships."))

rm(elo_combos, elo_combos_merged, elo_viz, elo_viz_res, con, o, temp_lm, a, i, n, var)

#---

#-West (WST)
#--Create dataframe of MMI scores and environmental data
wst_gam <- wst_full
wst_gam <-  wst_gam[(rownames(wst_gam) %in% rownames(chem_trans)), ]
wst_gam <- wst_gam[gtools::mixedorder(rownames(wst_gam)), ]
wst_global <- cbind(wst_gam, wst_chem)

#--Perform simple linear regressions for each environmental variable
wst_lm <- data.frame(Response = NA, Explanatory = NA, R2 = NA, p = NA)
for(i in 1:ncol(wst_chem)){
  a <- colnames(wst_chem)[i]
  temp_lm <- lm(wst_global$MMI ~ wst_global[[a]])
  wst_lm[i,1] <- "MMI"
  wst_lm[i,2] <- a
  wst_lm[i,3] <- summary(temp_lm)$adj.r.squared
  wst_lm[i,4] <- summary(temp_lm)$coefficients[2,4]
  r2 <- round(summary(temp_lm)$adj.r.squared, digits = 2)
  p <- round(summary(temp_lm)$coefficients[2,4], digits = 2)
  print(ggResidpanel::resid_panel(temp_lm))
  print(ggplot() +
          geom_point(aes(wst_global[[a]], wst_global$MMI), shape = 21, size = 2, color = "gray30", fill = "gray70", alpha = 0.50) +
          geom_smooth(aes(wst_global[[a]], wst_global$MMI), method = lm, se = TRUE, color = "gray30", fill = "gray", alpha = 0.50) +
          #labs(x = paste(a, "(transformed)"), y = "MMI score") +
          annotate(geom = "label", label = deparse(bquote(atop(paste(R^2 ~ "=" ~ .(r2)), italic(.("p") ~ "=" ~ .(p))))), parse = TRUE,
                   label.r = unit(0, "lines"), size = 3, x = Inf, y = Inf, hjust = 1, vjust = 1, label.size = 0.25, label.padding = unit(0.075, "lines")) +
          scale_x_continuous(expand = c(0, 0)) +
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank()))
  ggsave(paste("wst_lm_", a, ".png", sep = ""),
         plot = last_plot(),
         width = 2.4,
         height = 2.4,
         units = c("in"),
         dpi = 1200)
}
wst_lm

#--Perform multiple linear regression for all environmental variables
wst_mlm <- lm(MMI ~ ., data = wst_global[,8:15])
print(ggResidpanel::resid_panel(wst_mlm))
summary(wst_mlm)

#--Construct all possible GAMs with 1-2 explanatory variables
library(mgcv)

#---Create all 1-7 variable combinations using all variables
wst_combos <- list()
wst_combvec <- vector()

for(a in colnames(wst_chem)){
  value <- paste0("s(", a, ", k = 5)")
  wst_combvec[[a]] <- value
}
wst_combos$one <- as.matrix(t(combn(wst_combvec, 1)))
for(i in 1:nrow(wst_combos$one)){
  value <- paste(wst_combos$one[i,], collapse = " + ")
  wst_combos$one[i,1] <- value
}
wst_combos$one <- wst_combos$one[,-2]
wst_combos$two <- as.matrix(t(combn(wst_combvec, 2)))
for(i in 1:nrow(wst_combos$two)){
  value <- paste(wst_combos$two[i,], collapse = " + ")
  wst_combos$two[i,1] <- value
}
wst_combos$two <- wst_combos$two[,-2]
wst_combos$three <- as.matrix(t(combn(wst_combvec, 3)))
for(i in 1:nrow(wst_combos$three)){
  value <- paste(wst_combos$three[i,], collapse = " + ")
  wst_combos$three[i,1] <- value
}
wst_combos$three <- wst_combos$three[,-c(2:3)]
wst_combos$four <- as.matrix(t(combn(wst_combvec, 4)))
for(i in 1:nrow(wst_combos$four)){
  value <- paste(wst_combos$four[i,], collapse = " + ")
  wst_combos$four[i,1] <- value
}
wst_combos$four <- wst_combos$four[,-c(2:4)]
wst_combos$five <- as.matrix(t(combn(wst_combvec, 5)))
for(i in 1:nrow(wst_combos$five)){
  value <- paste(wst_combos$five[i,], collapse = " + ")
  wst_combos$five[i,1] <- value
}
wst_combos$five <- wst_combos$five[,-c(2:5)]
wst_combos$six <- as.matrix(t(combn(wst_combvec, 6)))
for(i in 1:nrow(wst_combos$six)){
  value <- paste(wst_combos$six[i,], collapse = " + ")
  wst_combos$six[i,1] <- value
}
wst_combos$six <- wst_combos$six[,-c(2:6)]
wst_combos$seven <- as.matrix(t(combn(wst_combvec, 7)))
for(i in 1:nrow(wst_combos$seven)){
  value <- paste(wst_combos$seven[i,], collapse = " + ")
  wst_combos$seven[i,1] <- value
}
wst_combos$seven <- wst_combos$seven[,-c(2:7)]
wst_combos_merged <- as.matrix(do.call(c, wst_combos))

rm(wst_combvec, value)

#---Add response variable name to variable combinations
wst_combos_merged <- as.matrix(paste("MMI ~", wst_combos_merged))

#---Run GAM all subsets
n <- nrow(wst_combos_merged)
wst_mods <- data.frame(formula = rep(NA, n),
                       max.pval = rep(NA, n),
                       dev.exp = rep(NA, n),
                       concurvity.worst = rep(NA, n),
                       concurvity.obs = rep(NA, n),
                       concurvity.est = rep(NA, n),
                       AIC = rep(NA, n))
for(i in 1:nrow(wst_combos_merged)){
  mod <- gam(as.formula(paste(wst_combos_merged[i,1])),
             data = wst_global,
             fx = FALSE,
             select = TRUE,
             gamma = 1.4,
             family = tw,
             method = "REML",
             na.action = "na.fail")
  con <- concurvity(mod)
  wst_mods[i, 1] <- paste0(summary(mod)$formula[2],summary(mod)$formula[1],summary(mod)$formula[3])
  wst_mods[i, 2] <- max(summary(mod)$s.pv)
  wst_mods[i, 3] <- summary(mod)$dev.expl
  wst_mods[i, 4] <- max(con[1, -1])
  wst_mods[i, 5] <- max(con[2, -1])
  wst_mods[i, 6] <- max(con[3, -1])
  wst_mods[i, 7] <- AIC(mod)
}

wst_stored <- wst_mods #Retain all subsets models for reference

#---Filter out nonsignificant, low deviance explained, and low performing GAMs
wst_mods <- wst_mods[wst_mods$max.pval <= 0.05, ]
wst_mods <- wst_mods[wst_mods$dev.exp >= 0.15, ]
wst_mods <- wst_mods[wst_mods$concurvity.est < 0.70, ]
wst_mods$concurvity.worst <- NULL
wst_mods$concurvity.obs <- NULL

#---Determine models with significant smooth terms or non-convergence in gam.check (also look at residual plots)
wst_mods$min_pval_check <- NA
wst_mods$convergence <- NA
for(i in 1:nrow(wst_mods)){
  mod <- gam(as.formula(paste(wst_mods[i,1])),
             data = wst_global,
             fx = FALSE, 
             select = TRUE,
             gamma = 1.4,
             family = tw,
             method = "REML",
             na.action = "na.fail")
  set.seed(10)
  wst_mods[i,6] <- min(k.check(mod)[,4]) #Minimum p-value of smooth terms in gam.check
  set.seed(10)
  wst_mods[i,7] <- capture.output(gam.check(mod))[3] #Whether model converged
}

#####Calculate delta AIC and Akaike weights and likelihoods
library(qpcR)
#wst_mods <- wst_mods[wst_mods$min_pval_check > 0.05, ] #Remove models with significant smooth terms in gam.check (i.e., oversmoothed). For MMI v24 results, removing these models results in no significant models. Therefore, this step is removed to allow for presentation of approximate results.
wst_mods <- wst_mods[grep("full convergence", wst_mods$convergence), ] #Remove re-run models that did not converge
wst_mods <- wst_mods[wst_mods$max.pval < 0.05, ]
wst_mods <- wst_mods[wst_mods$dev.exp >= 0.15, ]
wst_mods$delta.AIC <- akaike.weights(wst_mods$AIC)$deltaAIC
wst_mods$akaike.likelihoods <- akaike.weights(wst_mods$AIC)$rel.LL
wst_mods$akaike.weights <- akaike.weights(wst_mods$AIC)$weights
wst_mods <- wst_mods[wst_mods$delta.AIC < 2, ] #Retain only those models with delta AIC < 2
wst_mods <- wst_mods[order(-wst_mods$akaike.weights), ] #Order models by Akaike weight

#---Check GAM residuals to validate model
for(i in 1:nrow(wst_mods)){
  mod <- gam(as.formula(paste(wst_mods[i,1])),
             data = wst_global,
             fx = FALSE,
             select = TRUE,
             gamma = 1.4,
             family = tw,
             method = "REML",
             na.action = "na.fail")
  gam.check(mod)
}

#---Un-transform explanatory variable values for partial effects plots
df <- data.frame(transformation = c("quart", "cube", "square", "no transformation", "square root", "cube root", "fourth root", "log(x + 1)", "reciprocal root + 1", "reciprocal + 1", "reciprocal square +1", "log(x)", "reciprocal root", "reciprocal", "reciprocal square"),
                 formula = c("wst_viz[[i]][1]^0.25", "wst_viz[[i]][1]^(1/3)", "wst_viz[[i]][1]^0.5", "wst_viz[[i]][1]", "wst_viz[[i]][1]^2", "wst_viz[[i]][1]^3", "wst_viz[[i]][1]^4", "exp(wst_viz[[i]][1]) - 1", "(-1 / wst_viz[[i]][1])^2 - 1", "(-1 / wst_viz[[i]][1]) - 1", "(-1 / wst_viz[[i]][1])^0.5 - 1","exp(wst_viz[[i]][1])", "(-1 / wst_viz[[i]][1])^2", "-1 / wst_viz[[i]][1]", "(-1 / wst_viz[[i]][1])^0.5"),
                 blank = rep(NA, 15))
chem_trans_meta$wst_inverse <- lapply(chem_trans_meta$transformation, function(x) df$formula[match(x, df$transformation)])

#---Construct partial effects plots for filtered GAMs
#----Re-run model
mod <- gam(as.formula(paste(wst_mods[1,1])),
           data = wst_global,
           fx = FALSE,
           select = TRUE,
           gamma = 1.4,
           family = tw,
           method = "REML",
           na.action = "na.fail")

#----Extract fitted values and residuals
library(mgcViz)
wst_viz <- list()
wst_viz_res <- list()
for(i in 1:length(mod$var.summary)){
  var <- names(mod$var.summary)[[i]] 
  o <- plot(sm(getViz(mod), i))
  wst_viz[[var]] <- o$data$fit
  wst_viz_res[[var]] <- o$data$res
}

#----Process and plot fitted values and residuals
for(i in names(wst_viz)){
  a <- i
  wst_viz[[i]][[a]] <- wst_viz[[i]]$x
  wst_viz[[i]][[a]] <- eval(parse(text = paste0(as.data.frame(t(chem_trans_meta))[5,a])))
  wst_viz[[i]] <- as.data.frame(as.matrix(wst_viz[[i]]))
  
  wst_viz_res[[i]][[a]] <- wst_viz_res[[i]]$x
  wst_viz_res[[i]][[a]] <- eval(parse(text = gsub("viz", "viz_res", paste0(as.data.frame(t(chem_trans_meta))[5,a]))))
  wst_viz_res[[i]] <- as.data.frame(as.matrix(wst_viz_res[[i]]))
  
  plot(ggplot(wst_viz[[i]], aes(x = wst_viz[[i]][[a]], y = y)) +
         geom_ribbon(aes(ymin = y - se, ymax = y + se), fill = "gray", color = NA, alpha = 0.15) + #Plot confidence intervals
         #geom_point(data = wst_viz_res[[i]], aes(x = wst_viz_res[[i]][[a]], y = y), shape = 21, size = 2, color = "gray50", fill = "gray70") + #Optional: plot residuals that fall within fitted line
         geom_rug(data = wst_viz_res[[i]], aes(x = wst_viz_res[[i]][[i]], y = NULL), color = "gray30") + #Plot rugs along x-axis to represent distribution of all residuals
         geom_line(size = 1, color = "gray30") + #Plot fitted line
         #labs(title = paste("Model deviance explained = ", round(summary(mod)$dev.expl, digits = 2), sep = ""), x = a, y = "Partial effect") +
         scale_x_continuous(limits = c(min(wst_viz[[i]][[a]]), max(wst_viz[[i]][[a]])), expand = c(0, 0)) +
         scale_y_continuous(limits = c(min(wst_viz[[i]]$y - wst_viz[[i]]$se), max(wst_viz[[i]]$y + wst_viz[[i]]$se)), expand = c(0, 0)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank()))
  
  ggsave(paste("wst_gam_", i, ".png", sep = ""),
         plot = last_plot(),
         width = 2.4,
         height = 2.4,
         units = c("in"),
         dpi = 1200)
}

print(noquote("If 'wst_mods' has 0 observations, no GAMs are significant and explanatory. Disregard any plots or relationships."))

rm(wst_combos, wst_combos_merged, wst_viz, wst_viz_res, con, o, temp_lm, a, i, n, var)

#--------------

#Extract information for map making using ArcGIS
#-Extract latitude, longitude, and region and assign impairment based on CDPHE site MMI scores relative to reference site MMI scores used in the construction of the MMI
elo_map <- elo_full[,c(1,7)]
elo_map$Region <- "Plains and Lowlands"
elo_map$`MMI impairment` <- ifelse(elo_map$MMI < 53, "Likely Impacted", ifelse(elo_map$MMI > 66.5, "Likely Unimpacted", "Possibly Impacted")) #53 = national reference site minimum, 66.5 = national reference site 25th percentile

wst_map <- wst_full[,c(1,8)]
wst_map$Region <- "West"
wst_map$`MMI impairment` <- ifelse(wst_map$MMI < 37.5, "Likely Impacted", ifelse(wst_map$MMI > 55.5, "Likely Unimpacted", "Possibly Impacted")) #37.5 = national reference site minimum, 55.5 = national reference site 25th percentile

#-Merge elo_map and wst_map
map <- rbind(elo_map, wst_map)
map$`Region_biotic integrity` <- paste(map$Region, "_", map$`MMI impairment`, sep = "")
map$`Collection year` <- db[match(map$SAMPLEID, db$UID), "Sample Date"]
map$`Collection year` <- paste("FY", sub('.*(?=.{2}$)', '', map$`Collection year`, perl=T), sep = "")
map$`Collection date` <- db[match(map$SAMPLEID, db$UID), "Sample Date"]
map$`Sample type` <- db[match(map$SAMPLEID, db$UID), "Sample"]
map$`Station ID` <- db[match(map$SAMPLEID, db$UID), "StationID"]
map$`Waterbody name` <- db[match(map$SAMPLEID, db$UID), "WaterbodyName"]
map$`Collection location` <- db[match(map$SAMPLEID, db$UID), "Location"]
map$Latitude <- db[match(map$SAMPLEID, db$UID), "Lat_Dec"]
map$Longitude <- db[match(map$SAMPLEID, db$UID), "Long_Dec"]
map <- map[gtools::mixedorder(rownames(map)), ]
map <- merge(map, chem, by.x = "SAMPLEID", by.y = "UID", all.x = TRUE)
colnames(map)[1] <- "Sample ID"
colnames(map)[2] <- "MMI score"
colnames(map)[4] <- "MMI biotic integrity"
map[,13:19] <- NULL

#--Update environmental variable names
colnames(map)[13:33] <- c("pH", "Conductivity (uS/cm)", "DO", "Water temperature", "Alkalinity (mg/L)", "Dissolved arsenic (ug/L)",
                          "Dissolved copper (ug/L)", "Chloride (mg/L)", "Total hardness (mg/L)", "Dissolved iron (ug/L)",
                          "TREC iron (ug/L)", "Dissolved manganese (ug/L)", "Ammonia (mg/L)", "Nitrate + nitrite (mg/L)", "Kjeldahl nitrogen (mg/L)",
                          "Total nitrogen (mg/L)", "Total phosphorus (mg/L)", "Selenium (ug/L)", "Sulfate (mg/L)", "Zinc (ug/L)", 
                          "Total suspended solids (mg/L)")

#--Save map for upload to ArcGIS
#write.csv(map, "ArcGIS_Input_Map.csv", row.names = FALSE)

#--To create map:
#---Go to https://www.arcgis.com/home/webmap/viewer.html
#---Select Basemap > USGS National Map
#---Select Add > Browse Living Atlas Layers > Search "USA Counties" > Add ("+")
#----Adjust transparency, line color, etc. using Details > Content > USA Counties > Change Style
#---Select Add > Add Layer from File > "CDPHE_MMI_Map.csv"
#----To display MMI impairment by region, select Details > Choose an attribute to show: "Region impairment" > Select a drawing style: "Types (Unique symbols)" > adjust display options
#----To display other variables (e.g., phosphorus, conductivity), add new layer using same file ("CDPHE_MMI_Map.csv") and select target attribute

#--------------

#Perform qualitative precision analysis on MMI scores of duplicate samples
#-Calculate mean difference in MMI scores between each pair of samples
pairs <- as.data.frame(t(combn(nrow(map), 2)))
pairs$Difference <- NA
for(i in 1:nrow(pairs)){
  a <- pairs[i,1]
  b <- pairs[i,2]
  pairs[i,3] <- abs(map[a,2] - map[b,2])
}
pairs_agg <- aggregate(pairs, by = list(pairs$V1), FUN = "mean")
pairs_agg <- pairs_agg[,c(2,4)]
colnames(pairs_agg)[2] <- "Sample ID"
colnames(pairs_agg) <- c("Sample ID", "Mean difference")

#-Extract duplicate samples
map$Name <- paste(map$`Station ID`, map$`Collection date`, sep = "_")
dups <- map[duplicated(map$Name), ]
dups2 <- map[duplicated(map$Name, fromLast = TRUE), ]
dups <- rbind(dups, dups2)
dups <- dups[gtools::mixedorder(dups$Name), ]
dups$Difference <- NA
rm(dups2)

#-Calculate difference in MMI scores between duplicate samples
for(i in c(1,3,5,7,9,11,13,15,17,19)){
  dups[i,35] <- abs(dups[i,2] - dups[(i + 1),2])
}
for(i in c(2,4,6,8,10,12,14,16,18,20)){
  dups[i,35] <- abs(dups[i,2] - dups[(i - 1),2])
}

#-Create dataframe that has mean difference of all paired samples and difference of duplicate samples
pairs_agg$`Duplicate difference` <- dups[match(pairs_agg$`Sample ID`, dups$`Sample ID`), "Difference"]
pairs_agg <- pairs_agg[gtools::mixedorder(pairs_agg$`Duplicate difference`), ]
pairs_agg$`Station ID` <- dups[match(pairs_agg$`Sample ID`, dups$`Sample ID`), "Station ID"]
pairs_agg$`Collection Year` <- dups[match(pairs_agg$`Sample ID`, dups$`Sample ID`), "Collection year"]
pairs_agg <- pairs_agg[1:20,]

#--Save dataframe for appendix table
write.csv(pairs_agg, "Duplicate_Score_Differences.csv", row.names = FALSE)

#--------------

#Calculate summary statistics for MMI score distribution
mmi_stats <- data.frame(All = NA, West = NA, "Plains and Lowlands" = NA)
mmi_stats[1,3] <- nrow(elo_map)
mmi_stats[4,3] <- round(nrow(elo_map[grep("Likely Impacted", elo_map$`MMI impairment`),]) / nrow(elo_map) * 100, digits = 2)
mmi_stats[3,3] <- round(nrow(elo_map[grep("Possibly Impacted", elo_map$`MMI impairment`),]) / nrow(elo_map) * 100, digits = 2)
mmi_stats[2,3] <- round(nrow(elo_map[grep("Likely Unimpacted", elo_map$`MMI impairment`),]) / nrow(elo_map) * 100, digits = 2)

mmi_stats[1,2] <- nrow(wst_map)
mmi_stats[4,2] <- round(nrow(wst_map[grep("Likely Impacted", wst_map$`MMI impairment`),]) / nrow(wst_map) * 100, digits = 2)
mmi_stats[3,2] <- round(nrow(wst_map[grep("Possibly Impacted", wst_map$`MMI impairment`),]) / nrow(wst_map) * 100, digits = 2)
mmi_stats[2,2] <- round(nrow(wst_map[grep("Likely Unimpacted", wst_map$`MMI impairment`),]) / nrow(wst_map) * 100, digits = 2)

mmi_stats[1,1] <- nrow(map)
mmi_stats[4,1] <- round(nrow(map[grep("Likely Impacted", map$`MMI biotic integrity`),]) / nrow(map) * 100, digits = 2)
mmi_stats[3,1] <- round(nrow(map[grep("Possibly Impacted", map$`MMI biotic integrity`),]) / nrow(map) * 100, digits = 2)
mmi_stats[2,1] <- round(nrow(map[grep("Likely Unimpacted", map$`MMI biotic integrity`),]) / nrow(map) * 100, digits = 2)

rownames(mmi_stats) <- c("Samples", "Likely unimpacted", "Possibly impacted", "Likely impacted")
colnames(mmi_stats)[3] <- "Plains and Lowlands"

#--------------

#Create table of CDPHE taxonomy, MMI taxonomy, and whether taxon counts were used in MMI computation
lookup <- read.csv("CDPHE_Taxonomy_Lookup_Table.csv", header = TRUE, strip.white = TRUE)

cdphe_tax <- as.data.frame(cdphe_tax)
colnames(cdphe_tax) <- "CDPHE"
cdphe_tax$MMI <- lookup[match(cdphe_tax$CDPHE, lookup$CDPHE), "MMI"]
cdphe_tax$Revised <- ifelse(is.na(cdphe_tax$MMI), "N", "Y")
cdphe_tax$MMI <- ifelse(is.na(cdphe_tax$MMI), cdphe_tax$CDPHE, cdphe_tax$MMI)

#-Used in MMI?
missing_tax <- read.csv("MissingTaxaList.csv", header = TRUE, strip.white = TRUE)
cdphe_tax$`Used in MMI` <- missing_tax[match(cdphe_tax$MMI, missing_tax$TAXON), "TAXON"]
cdphe_tax$`Used in MMI` <- ifelse(is.na(cdphe_tax$`Used in MMI`), "Y", "N")

#-Fill in taxonomy lineage information
biodata <- read.csv("biodata_taxonomy_complete_20170614.csv", header = TRUE, strip.white = TRUE)
biodata$BiodataTaxonName <- gsub(" ", "_", biodata$BiodataTaxonName)
biodata$Form <- gsub("fo.", "f.", biodata$Form)
biodata <- biodata[!duplicated(biodata$BiodataTaxonName), ]
cdphe_tax <- merge(cdphe_tax, biodata[,c(10,22,25,29,33,36,38:41)], by.x = "MMI", by.y = "BiodataTaxonName", all.x = TRUE, all.y = FALSE)
cdphe_tax$Phylum <- gsub("Heterokontophyta", "Bacillariophyta", cdphe_tax$Phylum)
colnames(cdphe_tax)[7] <- "Order"

#-Save file
#write.csv(cdphe_tax, "MMI_CDPHE_Diatom_Taxonomy_Harmonization.csv", row.names = FALSE)
