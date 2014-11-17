###############################################################################
#			DECODE data
###############################################################################
# Kate Lee Feb 2014


### Overview
# ASSUMED that the default settings were used
# pvalue, beta and se scores were for the additive frequentist test (correspondence with BOWN)

# reduced data set with unix 
# (got RSID STRAND CHR POS EA NON_EA EAF_CASES BETA SE INFO)
#	cut -d ' ' -f 2,3,5,6,7,8,9,11,12,15 DECODE.AAA.All.noAdj.30JUL2013.GTh.txt > ../formatted/small_decode


# tidy in R
# 	re-name allele, and position columns (keep frequentist output names)
#	remove rows with indel data
#	removed rows with info <0.3
#	from correspondance with Maria:
#		removed rows with maf or eaf <0.01
#		We have beta and se
#		We need to
#		1.	Calculate OR= exp(beta)
#		2.	Calculate lower limit=beta-(1.96xse)
#		3.	Calculate upper limit= beta+(1.96xse)
#		4.	Calculate exp(lower limit)
#		5.	Calculate exp(upper limit)
#		In the final dataset we will include OR and what we calculated in 4 and 5.
# 	select columns required for analysis and write to new file

###   decodetidy.R file

# clean workspace
rm(list = ls(all=TRUE))


# read in snpfile 
snpfile <- read.table("../../formatted2/small_decode", header=TRUE)

# re-name some columns
colnames(snpfile)[colnames(snpfile) == "RSID"] <- "MARKER"
colnames(snpfile)[colnames(snpfile) == "EFFECT_ALLELE"] <- "EA"
colnames(snpfile)[colnames(snpfile) == "NON_EFFECT_ALLELE"] <- "NEA"
colnames(snpfile)[colnames(snpfile) == "EAF_CASES"] <- "EAF"
colnames(snpfile)[colnames(snpfile) == "INFO"] <- "info"


# change character class of useful columns 
snpfile[,"MARKER"] <- as.character(snpfile[,"MARKER"])
snpfile[,"CHR"] <- as.character(snpfile[,"CHR"])
snpfile[,"POS"] <- as.numeric(snpfile[,"POS"])
snpfile[,"EA"] <- as.character(snpfile[,"EA"])
snpfile[,"NEA"] <- as.character(snpfile[,"NEA"])
snpfile[,"EAF"] <- as.numeric(snpfile[,"EAF"])
#snpfile[,"PVAL"] <- as.numeric(snpfile[,"PVAL"])
snpfile[,"BETA"] <- as.numeric(snpfile[,"BETA"])
snpfile[,"SE"] <- as.numeric(snpfile[,"SE"])
snpfile[,"info"] <- as.numeric(snpfile[,"info"])

# remove rows with alleleA or alleleB that has a string longer than 1 character
snpfile <- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"EA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]
snpfile <- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"NEA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]

# remove rows where info is less than 0.3
snpfile <- snpfile[ snpfile$info > 0.3 , ]

# remove rows where maf or eaf is <0.01
snpfile <- snpfile[ snpfile$EAF > 0.01 , ]

# remove rows with missing data
snpfile <- na.omit(snpfile)

# rm rows where rsid = "."
snpfile <- snpfile[snpfile$MARKER != ".",]

# caluclate OR from BETA (OR=exp(beta))
snpfile$OR <- apply(snpfile[,"BETA", drop = FALSE], 1, exp)

# calculate log lower limit = beta - (1.96Xse)
snpfile[,"log_all_OR_lower"] <- snpfile[,"BETA"] - (snpfile[,"SE"]*1.96)

# calculate log upper limit = beta - (1.96Xse)
snpfile[,"log_all_OR_upper"] <- snpfile[,"BETA"] + (snpfile[,"SE"]*1.96)

# calculate lower limit = exp(log_all_OR_lower)
snpfile$OR_95L <- apply(snpfile[,"log_all_OR_lower", drop = FALSE], 1, exp)

# calculate upper limit = exp(log_all_OR_upper)
snpfile$OR_95U <- apply(snpfile[,"log_all_OR_upper", drop = FALSE], 1, exp)

# order by rsid
snpfile <- snpfile[order(snpfile$MARKER),]

# subset data
newsnpfile <- snpfile[,c("MARKER", "CHR", "POS", "NEA", "EA", "OR", "OR_95L", "OR_95U", "EAF")]

#write to new csv file
write.csv(newsnpfile, "../../formatted2/DECODE.AAA.All.noAdj.30JUL2013.GTh.csv", row.names=FALSE)

