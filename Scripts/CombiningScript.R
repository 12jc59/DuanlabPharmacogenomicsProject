
library(readr)
European_LD_sorted <- read_delim("European_LD_sorted.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
EastAsian_LD_sorted <- read_delim("EastAsian_LD_sorted.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
American_LD_sorted <- read_delim("American_LD_sorted.txt", 
                                 "\t", escape_double = FALSE, col_names = FALSE, 
                                 trim_ws = TRUE)
African_LD_sorted <- read_delim("African_LD_sorted.txt", 
                                  "\t", escape_double = FALSE, col_names = FALSE, 
                                  trim_ws = TRUE)



EastAsianPairs = c()
for (i in 1:dim(EastAsian_LD_sorted)[1]){
  EastAsianPairs[i] <- paste(EastAsian_LD_sorted$X1[i],EastAsian_LD_sorted$X2[i])
}

EuropeanPairs = c()
for (i in 1:dim(European_LD_sorted)[1]){
  EuropeanPairs[i] <- paste(European_LD_sorted$X1[i],European_LD_sorted$X2[i])
}

AmericanPairs = c()
for (i in 1:dim(American_LD_sorted)[1]){
  AmericanPairs[i] <- paste(American_LD_sorted$X1[i],American_LD_sorted$X2[i])
}

AfricanPairs = c()
for (i in 1:dim(African_LD_sorted)[1]){
  AfricanPairs[i] <- paste(African_LD_sorted$X1[i],African_LD_sorted$X2[i])
}


CombinedPairs<-c(EastAsianPairs,EuropeanPairs,AfricanPairs,AmericanPairs)
CombinedPairs<-unique(CombinedPairs)

FINAL_TABLE<-data.frame(MyVariant=character(),#1
                        MyVariant_chromosome=character(),#2
                        MyVariant_gene=character(),#3
                        
                        MyVariant_biotype=character(),#4
                        MyVariant_Impact=character(),#5
                        MyVariant_Annotation=character(),#6
                        MyVariant_RegulomeDB_score=character(),#7
                        
                        MyVariant_SNPedia_Suumary=character(),#8
                        MyVariant_SNPedia_Trait=character(),#9
                        MyVariant_SNPedia_Disease=character(),#10
                        MyVariant_SNPedia_Drugclasses=character(),#11
                        MyVariant_SNPedia_Annotation=character(),#12
                        
                        
                        PharmGKBVariant=character(),#13
                        PharmGKBVariant_chromosome=character(),#14
                        PharmGKBVariant_gene=character(),#15
                        
                        PharmGKBVariant_biotype=character(),#16
                        PharmGKBVariant_Impact=character(),#17
                        PharmGKBVariant_Annotation=character(),#18
                        PharmGKBVariant_RegulomeDB_score=character(),#19
                        
                        PharmGKBVariant_SNPedia_Suumary=character(),#20
                        PharmGKBVariant_SNPedia_Trait=character(),#21
                        PharmGKBVariant_SNPedia_Disease=character(),#22
                        PharmGKBVariant_SNPedia_Drugclasses=character(),#23
                        PharmGKBVariant_SNPedia_Annotation=character(),#24
                        
                        European_r2=integer(),#25
                        EastAsian_r2=integer(),#26
                        American_r2=integer(),#27
                        African_r2=integer(),#28
                        stringsAsFactors=FALSE)
# writing rsID & LD information
for (i in 1:length(CombinedPairs)){
  # rsid
  FINAL_TABLE[i,1]<-strsplit(CombinedPairs[[i]]," ")[[1]][1]
  FINAL_TABLE[i,13]<-strsplit(CombinedPairs[[i]]," ")[[1]][2]
  
  if (is.na(match(CombinedPairs[i],EuropeanPairs))==FALSE){
    FINAL_TABLE[i,25]<-European_LD_sorted$X3[match(CombinedPairs[i],EuropeanPairs)]}
  
  if (is.na(match(CombinedPairs[i],EastAsianPairs))==FALSE){
    FINAL_TABLE[i,26]<-EastAsian_LD_sorted$X3[match(CombinedPairs[i],EastAsianPairs)]}
  
  if (is.na(match(CombinedPairs[i],AmericanPairs))==FALSE){
    FINAL_TABLE[i,27]<-American_LD_sorted$X3[match(CombinedPairs[i],AmericanPairs)]}
  
  if (is.na(match(CombinedPairs[i],AfricanPairs))==FALSE){
    FINAL_TABLE[i,28]<-African_LD_sorted$X3[match(CombinedPairs[i],AfricanPairs)]}
}

# read in snpeff output
OUT_combined_cleaned <- read_csv("C:/Users/J/Desktop/pharmGpaper/snpEff/OUT_combined_cleaned.csv")

for (i in 1:length(CombinedPairs)){
  #for Myvariants
  FINAL_TABLE[i,2]<-OUT_combined_cleaned$CHR[match(FINAL_TABLE[i,1],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,3]<-OUT_combined_cleaned$Gene[match(FINAL_TABLE[i,1],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,4]<-OUT_combined_cleaned$Biotype[match(FINAL_TABLE[i,1],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,5]<-OUT_combined_cleaned$Impact[match(FINAL_TABLE[i,1],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,6]<-OUT_combined_cleaned$Prediction[match(FINAL_TABLE[i,1],OUT_combined_cleaned$RS)]

  #for PharmGKBvariants
  FINAL_TABLE[i,14]<-OUT_combined_cleaned$CHR[match(FINAL_TABLE[i,13],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,15]<-OUT_combined_cleaned$Gene[match(FINAL_TABLE[i,13],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,16]<-OUT_combined_cleaned$Biotype[match(FINAL_TABLE[i,13],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,17]<-OUT_combined_cleaned$Impact[match(FINAL_TABLE[i,13],OUT_combined_cleaned$RS)]
  FINAL_TABLE[i,18]<-OUT_combined_cleaned$Prediction[match(FINAL_TABLE[i,13],OUT_combined_cleaned$RS)]
}


#Load regulomeDB data
regulomDB_cleaned <- read_csv("regulomDB_results.csv", 
                              col_names = FALSE)

for (i in 1:length(CombinedPairs)){
  #for Myvariants
  FINAL_TABLE[i,7]<-regulomDB_cleaned$X2[match(FINAL_TABLE[i,1],regulomDB_cleaned$X1)]
  FINAL_TABLE[i,19]<-regulomDB_cleaned$X2[match(FINAL_TABLE[i,13],regulomDB_cleaned$X1)]
}




### SNPedia search
# to only grab rs# that exists on SNPedia - if unfound SNP is throwin in parse function, it crashes
snpset<-getPages (titles = OUT_combined_cleaned$RS)
snpset<-Filter(Negate(is.null), snpset)

pg<-getPages (
  titles = names(snpset),
  wikiParseFunction = extractGenotypeTags,
  tags = c ("Summary","Trait","Diseases","Drug Classes","Annotation")
)
variant_summary<-as.data.frame(pg)

#Write out and reformat the summary table
#write_tsv(variant_summary,"OUT")

variants_with_phenotype <- read_delim("variants_with_phenotype.txt", 
                                      "\t", escape_double = FALSE, trim_ws = TRUE)

for (i in 1:length(CombinedPairs)){
  #for Myvariants
  FINAL_TABLE[i,8]<-variants_with_phenotype$Summary[match(FINAL_TABLE[i,1],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,9]<-variants_with_phenotype$Trait[match(FINAL_TABLE[i,1],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,10]<-variants_with_phenotype$Diseases[match(FINAL_TABLE[i,1],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,11]<-variants_with_phenotype$Drug_Classes[match(FINAL_TABLE[i,1],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,12]<-variants_with_phenotype$Annotation[match(FINAL_TABLE[i,1],variants_with_phenotype$Variant)]
  
  #for PharmGKBvariants
  FINAL_TABLE[i,20]<-variants_with_phenotype$Summary[match(FINAL_TABLE[i,13],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,21]<-variants_with_phenotype$Trait[match(FINAL_TABLE[i,13],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,22]<-variants_with_phenotype$Diseases[match(FINAL_TABLE[i,13],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,23]<-variants_with_phenotype$Drug_Classes[match(FINAL_TABLE[i,13],variants_with_phenotype$Variant)]
  FINAL_TABLE[i,24]<-variants_with_phenotype$Annotation[match(FINAL_TABLE[i,13],variants_with_phenotype$Variant)]
}


write_tsv(FINAL_TABLE,"FINAL_TABLE")