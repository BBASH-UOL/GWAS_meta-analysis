###############################################################################
#                       BBASH DATA  
###############################################################################


### Overview
# was originally run in snptest by chromosome
# concatenated these using mergebbashdata.PBS bash script
# 	this took out the comment and header lines in all files but chr1
#	added in the chromosome information
#	concatenated each file in order to the chr1 file
#	resulted in a file called snptest_cov3_all_chr with all the data

# used unix to pare down file 
# (got rsid chromosome position alleleA alleleB info cases_maf all_OR all_OR_lower all_OR_upper frequentist_add_pvalue frequentist_add_info frequentist_add_beta_1 frequentist_add_se_1)
# note, kept frequentist beta and se to be able to compare it with all_OR

#	tail -n +13 snptest_cov3_all_chr > headless_snptest_cov3_all_chr
#	cut -d ' ' -f 2,3,4,5,6,9,30,39,40,41,42,43,44,45 headless_snptest_cov3_all_chr > ../formatted/small_snptest
#	rm headless_snptest_cov3_all_chr
	

# tidying in R (tested on first 200 lines of file)
#	renamed columns to GWAMA-friendly names
#	removed rows with indel data
# 	add re-calculated beta and se columns (logOR, se_logOR) to check against beta_1 and se_1
#	renamed alleleB as effect allele (confirmed by Maria with Matt Bown)
#	recalcuated adj se and upper and lower bounds for chr 7,14,15,16,X
#	removed rows with info <0.3
#	removed rows with maf or eaf <0.01
# 	select columns required for analysis and write to new file
#	select columns required to check data and write to a second file


###    bbashtidy.R file

# clean workspace
rm(list = ls(all=TRUE))

# read in snpfile
snpfile <- read.table("../../formatted2/small_snptest", header=TRUE)

# re-name some columns
colnames(snpfile)[colnames(snpfile) == "alleleA"] <- "NEA"
colnames(snpfile)[colnames(snpfile) == "alleleB"] <- "EA"
colnames(snpfile)[colnames(snpfile) == "rsid"] <- "MARKER"
colnames(snpfile)[colnames(snpfile) == "chromosome"] <- "CHR"
colnames(snpfile)[colnames(snpfile) == "position"] <- "POS"
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

# remove rows with alleleA or alleleB that has a string longer than 1 character
snpfile<- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"NEA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]
snpfile <- snpfile[grep("^A$|^G$|^T$|^C$", snpfile[,"EA"], ignore.case = TRUE, perl = TRUE, value = FALSE),]

# remove rows where info is less than 0.3
snpfile <- snpfile[ snpfile$info > 0.3 , ]

# remove rows where maf or eaf is <0.01
snpfile <- snpfile[ snpfile$EAF > 0.01 , ]


# adjust se for chr 7,14,15,16,X (numbers from Maria Viskaduraki)
#	Adj se = se x sqrt(lambda)
#	Chromosome 7, lambda=1.113
#	Chrom 14, lambda=1.127
#	Chrom 15, lambda=1.124
#	Chrom 16, lambda=1.118
#	Chrom X, lambda=2.02

snpfile$adjse <- ifelse((snpfile$CHR == 7), snpfile$frequentist_add_se_1 * 1.054988, snpfile$frequentist_add_se_1)
snpfile$adjse <- ifelse((snpfile$CHR == 14), snpfile$adjse * sqrt(1.127), snpfile$adjse) 
snpfile$adjse <- ifelse((snpfile$CHR == 15), snpfile$adjse * sqrt(1.124), snpfile$adjse) 
snpfile$adjse <- ifelse((snpfile$CHR == 16), snpfile$adjse * sqrt(1.118), snpfile$adjse) 
snpfile$adjse <- ifelse((snpfile$CHR == "X"), snpfile$adjse * sqrt(2.02), snpfile$adjse) 

### adjust upr and lower limit for chr 7,14,15,16,X
#Then lower limit=OR-(1.96 x adj se)
snpfile$OR_95L <- ifelse((snpfile$CHR == 7 | snpfile$CHR  == 14 | snpfile$CHR  == 15 | snpfile$CHR  == 16 | snpfile$CHR  == "X"), snpfile$OR - (snpfile$adjse * 1.96), snpfile$OR_95L) 

#And upper limit=OR +(1.96 x adj se)
snpfile$OR_95U <- ifelse((snpfile$CHR  == 7 | snpfile$CHR  == 14 | snpfile$CHR  == 15 | snpfile$CHR  == 16 | snpfile$CHR  == "X"), snpfile$OR + (snpfile$adjse * 1.96), snpfile$OR_95U)


# order by rsid
snpfile <- snpfile[order(snpfile$MARKER),]

# subset data
newsnpfile <- snpfile[,c("MARKER", "CHR", "POS", "NEA", "EA", "OR", "OR_95L", "OR_95U", "EAF")]

# remove rows with missing data
newsnpfile <- na.omit(newsnpfile)

#write to new csv file
write.csv(newsnpfile, "../../formatted2/snptest_cov3_all_chr_formatted.csv", row.names=FALSE)
