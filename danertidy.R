
### Overview
# ASSUMED that the default settings were used
## BOWN correspondance:

#For the Dutch dataset the output is standard PLINK dosage output:

#     CHR   Chromosome code, if map file specified
#     SNP   SNP code
#      BP   Base-pair position, if map file specified
#      A1   Allele 1 code
#      A2   Allele 2 code
#     FRQ   Frequency of A1, from dosage data
#    INFO   R-squared quality metric / information content
#      OR   Odds ratio for association (or BETA for quantitative traits)
#      SE   Standard error of effect estimate
#       P   p-value for association tests

#Therefore the OR is the OR and SE is SE OR.


# get subset of file from unix
# got  CHR SNP BP  A1  A2   FRQ_A_813 INFO OR SE
#	cat daner_ane_aaa1_eur.ch.fl | tr -s ' '| cut -d ' ' -f 2,3,4,5,6,7,9,10,11 > ../formatted/small_daner


# tidy in R
# 	re-name allele, and position columns (keep frequentist output names)
#	remove rows with indel data
#	removed rows with info <0.3
#	removed rows with maf or eaf <0.01
#	from correspondance with Maria:
#		We have OR and se.
#		We need to 
#		1.	convert OR to logOR. 
#		2.	Calculate se(logOR)=se/OR
#		3.	Calculate lower limit=logOR-[1.96xse(logOR)]
#		4.	Calculate upper limit= logOR+[1.96xse(logOR)]
#		5.	Calculate exp(lower limit)
#		6.	Calculate exp(upper limit)
#		In the final dataset we will include OR and what we calculated in 5 and 6.
# 	select columns required for analysis and write to new file

###   danertidy.R file


# clean workspace
rm(list = ls(all=TRUE))


# read in snpfile
snpfile <- read.table("../../formatted2/small_daner", header=TRUE)

# re-name some columns
colnames(snpfile)[colnames(snpfile) == "SNP"] <- "MARKER"
colnames(snpfile)[colnames(snpfile) == "A1"] <- "EA"
colnames(snpfile)[colnames(snpfile) == "A2"] <- "NEA"
colnames(snpfile)[colnames(snpfile) == "INFO"] <- "info"
colnames(snpfile)[colnames(snpfile) == "BP"] <- "POS"
colnames(snpfile)[colnames(snpfile) == "FRQ_A_813"] <- "EAF"  ####### not sure about this??????

# change character class of useful columns 
snpfile[,"MARKER"] <- as.character(snpfile[,"MARKER"])
snpfile[,"CHR"] <- as.character(snpfile[,"CHR"])
snpfile[,"POS"] <- as.numeric(snpfile[,"POS"])
snpfile[,"EA"] <- as.character(snpfile[,"EA"])
snpfile[,"NEA"] <- as.character(snpfile[,"NEA"])
snpfile[,"EAF"] <- as.numeric(snpfile[,"EAF"])
snpfile[,"OR"] <- as.numeric(snpfile[,"OR"])
#snpfile[,"PVAL"] <- as.numeric(snpfile[,"P"])
snpfile[,"SE"] <- as.numeric(snpfile[,"SE"])

# remove rows with alleleA or alleleB that has a string longer than 1 character 
snpfile <- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"EA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]
snpfile <- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"NEA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]

# remove rows where info is less than 0.3
snpfile <- snpfile[ snpfile$info > 0.3 , ]

# remove rows where maf or eaf is <0.01
snpfile <- snpfile[ snpfile$EAF > 0.01 , ]

# remove rows with missing data
snpfile <- na.omit(snpfile)

# re-calculate logOR 
snpfile[,"logOR"] <- apply(snpfile[,"OR", drop = FALSE], 1, log)

# calculate se(logOR) = se/OR
snpfile$se_logOR <- snpfile[,"SE"]/snpfile[,"OR"]

# calculate log lower limit = logOR-(1.96*se_logOR)
snpfile$log_all_OR_lower <- snpfile[,"logOR"] - (snpfile[,"se_logOR"]*1.96)

# calculate log upper limit = logOR+(1.96*se_logOR)
snpfile$log_all_OR_upper <- snpfile[,"logOR"] + (snpfile[,"se_logOR"]*1.96)

# calculate lower limit = exp(log_all_OR_lower)
snpfile$OR_95L <- apply(snpfile[,"log_all_OR_lower", drop = FALSE], 1, exp)

# calculate upper limit = exp(log_all_OR_upper)
snpfile$OR_95U <- apply(snpfile[,"log_all_OR_upper", drop = FALSE], 1, exp)

# order by rsid
snpfile <- snpfile[order(snpfile$MARKER),]

# subset data
newsnpfile <- snpfile[,c("MARKER", "CHR", "POS", "NEA", "EA", "OR", "OR_95L", "OR_95U", "EAF")]

#write to new csv file
write.csv(newsnpfile, "../../formatted2/daner_ane_aaa1_eur.ch.csv", row.names=FALSE)


