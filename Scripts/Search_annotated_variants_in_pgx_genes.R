
PharmGKB_annotated_variants <- read.csv("PharmGKB_annotated_rsid.tsv",sep = "\t",as.is = TRUE)
Pgx160_Genomic_region <- read.csv("160Pgx_genomic_regions.csv",as.is = TRUE)

list=data.frame()
counter = 1
for (variant in 1:dim(PharmGKB_annotated_variants)[1]){
  for (geneLocation in 1:dim(Pgx160_Genomic_region)[1]){
    if (PharmGKB_annotated_variants[variant,]$Chromosome == Pgx160_Genomic_region[geneLocation,]$chrom){
      if(PharmGKB_annotated_variants[variant,]$Position_GRCh37 > Pgx160_Genomic_region[geneLocation,]$Start_padded.5000){
        if(PharmGKB_annotated_variants[variant,]$Position_GRCh37 < Pgx160_Genomic_region[geneLocation,]$End_padded.5000){
          list[counter,1]<-PharmGKB_annotated_variants[variant,]$RSID
          list[counter,2]<-PharmGKB_annotated_variants[variant,]$Gene.Symbols
          list[counter,3]<-Pgx160_Genomic_region[geneLocation,]$Name2
          counter = counter +1
        }
      }
    }
  }
}

write.csv(list,"PharmGKB_within160pgx.csv")
