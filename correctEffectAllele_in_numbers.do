*this file checks if the same effect allele was used for each snp*
*it also takes the appropriate action to correct alleles accordingly*



gen logOR_y=ln(OR_y)
gen loglb=ln(OR_95L_y)
gen logub=ln(OR_95U_y)
gen change=5




**************for Danner ONLY***
gen newNEA_y=4 if NEA_y==3
replace newNEA_y=3 if NEA_y==2
replace newNEA_y=2 if NEA_y==1

label define newNEA_y 1 "A" 2 "C" 3 "G" 4 "T"
label values newNEA_y newNEA_y
drop NEA_y
rename newNEA_y NEA_y

gen newEA_y=EA_yDanner
replace newEA_y=4 if EA_y==3
label define newEA_y 1 "A" 2 "C" 3 "G" 4 "T"
label values newEA_y newEA_y
drop EA_y
rename newEA_y EA_y



*for NON-complementary bases GWAMA does the work!!!! I have not used the code that follows*
**Go to complementary bases!!!******

replace change=0 if NEA_x==NEA_y & EA_x==EA_y
*-------------------------------------------*
***I have replaced A=1 C=2 G=3 T=4*****

replace change=2 if NEA_x==1 & EA_x==3 & NEA_y==3 & EA_y==1
replace change=2 if NEA_x==3 & EA_x==1 & NEA_y==1 & EA_y==3

**ACTIONch1**
replace change=1 if NEA_x==1 & EA_x==3 & NEA_y==4 & EA_y==2 
replace change=1 if NEA_x==3 & EA_x==1 & NEA_y==2 & EA_y==4 

replace change=2 if NEA_x==1 & EA_x==3 & NEA_y==2 & EA_y==4
replace change=2 if NEA_x==2 & EA_x==1 & NEA_y==4 & EA_y==2
*----------------------------------*
replace logORnew_y=-logOR_y if NEA_x==1 & EA_x==2 & NEA_y==2 & EA_y==1
replace logORnew_y=-logOR_y if NEA_x==2 & EA_x==1 & NEA_y==1 & EA_y==2
**NO ACTION**
replace change=1 if NEA_x==1 & EA_x==2 & NEA_y==4 & EA_y==3 
replace change=1 if NEA_x==2 & EA_x==1 & NEA_y==3 & EA_y==4 

replace logORnew_y=-logOR_y if NEA_x==1 & EA_x==2 & NEA_y==3 & EA_y==4
replace logORnew_y=-logOR_y if NEA_x==2 & EA_x==1 & NEA_y==4 & EA_y==3
*-------------------------------------------*
replace logORnew_y=-logOR_y if NEA_x==4 & EA_x==2 & NEA_y==2 & EA_y==4
replace logORnew_y=-logOR_y if NEA_x==2 & EA_x==4 & NEA_y==4 & EA_y==2
**NO ACTION**
replace change=1 if NEA_x==4 & EA_x==2 & NEA_y==1 & EA_y==3 
replace change=1 if NEA_x==2 & EA_x==4 & NEA_y==3 & EA_y==1 

replace logORnew_y=-logOR_y if NEA_x==4 & EA_x==2 & NEA_y==3 & EA_y==1
replace logORnew_y=-logOR_y if NEA_x==2 & EA_x==4 & NEA_y==1 & EA_y==3
*----------------------------------------------*
replace logORnew_y=-logOR_y if NEA_x==4 & EA_x==3 & NEA_y==3 & EA_y==4
replace logORnew_y=-logOR_y if NEA_x==3 & EA_x==4 & NEA_y==4 & EA_y==3
**NO ACTION**
replace change=1 if NEA_x==4 & EA_x==3 & NEA_y==1 & EA_y==2 
replace change=1 if NEA_x==3 & EA_x==4 & NEA_y==2 & EA_y==1 

replace change=2 if NEA_x==4 & EA_x==3 & NEA_y==2 & EA_y==1
replace change=2 if NEA_x==3 & EA_x==4 & NEA_y==1 & EA_y==2


**-----------------------------------------------------------------------------------------------------------------***
*for complementary bases*
** AT***
replace change=1 if NEA_x==1 & EA_x==4 & NEA_y==4 & EA_y==1 & EAF_x>0.55 & EAF_y>0.55
replace change=1 if NEA_x==4 & EA_x==1 & NEA_y==1 & EA_y==4 & EAF_x>0.55 & EAF_y>0.55
replace change=1 if NEA_x==1 & EA_x==4 & NEA_y==4 & EA_y==1 & EAF_x<0.45 & EAF_y<0.45
replace change=1 if NEA_x==4 & EA_x==1 & NEA_y==1 & EA_y==4 & EAF_x<0.45 & EAF_y<0.45

replace change=2 if NEA_x==1 & EA_x==4 & NEA_y==4 & EA_y==1 & EAF_x>0.55 & EAF_y<0.45
replace change=2 if NEA_x==4 & EA_x==1 & NEA_y==1 & EA_y==4 & EAF_x>0.55 & EAF_y<0.45
replace change=2 if NEA_x==1 & EA_x==4 & NEA_y==4 & EA_y==1 & EAF_x<0.45 & EAF_y>0.55
replace change=2 if NEA_x==4 & EA_x==1 & NEA_y==1 & EA_y==4 & EAF_x<0.45 & EAF_y>0.55

***CG now***Danner
replace change=1 if NEA_x==2 & EA_x==3 & NEA_y==3 & EA_y==2 & EAF_x>0.55 & EAF_y>0.55
replace change=1 if NEA_x==3 & EA_x==2 & NEA_y==2 & EA_y==3 & EAF_x>0.55 & EAF_y>0.55
replace change=1 if NEA_x==2 & EA_x==3 & NEA_y==3 & EA_y==2 & EAF_x<0.45 & EAF_y<0.45
replace change=1 if NEA_x==3 & EA_x==2 & NEA_y==2 & EA_y==3 & EAF_x<0.45 & EAF_y<0.45

replace change=2 if NEA_x==2 & EA_x==3 & NEA_y==3 & EA_y==2 & EAF_x>0.55 & EAF_y<0.45
replace change=2 if NEA_x==3 & EA_x==2 & NEA_y==2 & EA_y==3 & EAF_x>0.55 & EAF_y<0.45
replace change=2 if NEA_x==2 & EA_x==3 & NEA_y==3 & EA_y==2 & EAF_x<0.45 & EAF_y>0.55
replace change=2 if NEA_x==3 & EA_x==2 & NEA_y==2 & EA_y==3 & EAF_x<0.45 & EAF_y>0.55



***changes that I make if change==1 NO ACTION other than switching names***
replace NEA_y=NEA_x if change==1
replace EA_y=EA_x if change==1
***changes I make if change==2, that is alleles switched OR alleles switched and wrong strand***Untitled5
replace NEA_y=NEA_x if change==2
replace EA_y=EA_x if change==2

replace EAF_y=1-EAF_y if change==2
replace OR_y=1/OR_y if change==2
replace logOR_y=-logOR_y if change==2
replace OR_95L_y=exp(-loglb) if change==2
replace OR_95U_y=exp(-logub) if change==2

drop if EAF_y<0.55 & EAF_y>0.45 & change==1
drop if EAF_y<0.55 & EAF_y>0.45 & change==2

drop  NEA_x EA_x OR_x OR_95L_x OR_95U_x EAF_x CHR_x POS_x logOR_y loglb logub change

rename CHR_y CHR
rename POS_y POS
rename OR_y OR
rename OR_95L_y OR_95L
rename OR_95U_y OR_95U
rename EAF_y EAF
rename NEA_y NEA
rename EA_y EA

*******************************************************************************************

twoway (scatter OR_y MARKER) (spike OR_95L_y MARKER) (spike OR_95U_y MARKER) if OR_y<10
twoway (scatter OR_y MARKER) (spike OR_95L_y MARKER) (spike OR_95U_y MARKER) if OR_y>10

