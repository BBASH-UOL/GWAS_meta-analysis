###############################################################################
#			Geisinger data
##############################################################################
# Kate Lee 2014

### Overview
# seem to have used defaults in snp test, but have a few different naming schemes
# compared with bbash set
# merged chromosome files Matt created in order to input chromosome data using mergeGS.PBS

# used unix to pare down file 
# (got MARKER chromosome position alleleA alleleB cases_maf all_OR all_OR_lower all_OR_upper frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1)
# note, kept frequentist beta and se to be able to compare it with all_OR
#	cut -d , -f 2,3,4,5,6,9,25,34,35,36,37,38,39,40 GeisingerSNPtestResult_CHR.txt > ../formatted/small_geisinger


# tidying in R (tested on first 200 lines of file)
# 	re-name allele, and position columns (keep frequentist output names)
# 	add re-calculated logOR and se_logOR
#	remove rows with indel data
#	removed rows with info <0.3
#	removed rows with maf or eaf <0.01
# 	select columns required for analysis and write to new file
#	select columns required to check data and write to a second file

###   gstidy.R file

# clean workspace
rm(list = ls(all=TRUE))

# read in snpfile 
snpfile <- read.table("../../formatted2/small_geisinger", header=TRUE)


# re-name some columns
colnames(snpfile)[colnames(snpfile) == "allele_A"] <- "NEA"
colnames(snpfile)[colnames(snpfile) == "allele_B"] <- "EA"
colnames(snpfile)[colnames(snpfile) == "pos"] <- "POS"
colnames(snpfile)[colnames(snpfile) == "rsid"] <- "MARKER"
colnames(snpfile)[colnames(snpfile) == "chromosome"] <- "CHR"
colnames(snpfile)[colnames(snpfile) == "all_OR"] <- "OR"
colnames(snpfile)[colnames(snpfile) == "all_OR_lower"] <- "OR_95L"
colnames(snpfile)[colnames(snpfile) == "all_OR_upper"] <- "OR_95U"
colnames(snpfile)[colnames(snpfile) == "cases_maf"] <- "EAF"

# change character class of useful columns 
snpfile[,"MARKER"] <- as.character(snpfile[,"MARKER"])
snpfile[,"CHR"] <- as.character(snpfile[,"CHR"])
snpfile[,"POS"] <- as.numeric(snpfile[,"POS"])
snpfile[,"NEA"] <- as.character(snpfile[,"NEA"])
snpfile[,"EA"] <- as.character(snpfile[,"EA"])
snpfile[,"info"] <- as.character(snpfile[,"info"])
snpfile[,"EAF"] <- as.numeric(snpfile[,"EAF"])
snpfile[,"OR"] <- as.numeric(snpfile[,"OR"])
snpfile[,"OR_95L"] <- as.numeric(snpfile[,"OR_95L"])
snpfile[,"OR_95U"] <- as.numeric(snpfile[,"OR_95U"])

# remove rows with NEA or EA that has a string longer than 1 character (indel info)
newtmp <- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"NEA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]
newsnpfile <- newtmp[grep("^A$|^G$|^T$|^C$", newtmp[,"EA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]

# remove rows where info is less than 0.3
snpfile <- snpfile[ snpfile$info > 0.3 , ]

# remove rows where maf or eaf is <0.01
snpfile <- snpfile[ snpfile$EAF > 0.01 , ]

# remove rows with missing data
snpfile <- na.omit(snpfile)

# order by rsid
snpfile <- snpfile[order(snpfile$MARKER),]

# subset data
newsnpfile <- snpfile[,c("MARKER", "CHR", "POS", "NEA", "EA", "OR", "OR_95L", "OR_95U", "EAF")]

#write to new csv file
write.csv(newsnpfile, "test", row.names=FALSE)


