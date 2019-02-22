library(dplyr)
library(readr)
library(magrittr)
library(data.table)
library(reshape2)

#gyrA and parC mutations
gyrA_FQR <- c('S83L', 'D87N', 'D87L')
parC_FQR <- c('S80I', 'S80R', 'E84G', 'E84K', 'E84V')

#Detection of mutations in the gyrA and parC
#genes in Escherichia coli isolates carrying
#plasmid-mediated quinolone resistance genes
#from diseased food-producing animals
#There were six types of mutation
#found in parC: 132 (89.8 % of the 147 mutated isolates)
#isolates had S80I, five (3.4 %) had S80R, six (4.1 %) had
#E84K, two (1.4 %) had S80I and E84K, one (0.7 %) had
#S80I and E84G, and one (0.7 %) had S80I and A108V.

#mutations <- "gyrA: Ser-83-Leu, Asp-87-Asn; parC: Ser-80-Ile, Glu-84-Gly"
#https://github.com/katholt/Kleborate#resistance-gene-detection

#read in ariba output
gyrA_parC <- read_delim("Liv_gyrA_parC.tsv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

#change column names
colnames(gyrA_parC) <- gsub(pattern = "#", replacement = "", colnames(gyrA_parC))
colnames(gyrA_parC) <- c(colnames(gyrA_parC)[1:(ncol(gyrA_parC)-1)], "name")

#generate dataframe consisting of all samples in collection
name_list <- unique(gyrA_parC$name)
name_list <- name_list %>% as.data.frame()
colnames(name_list) <- 'name'

#remove repeated colnames and change column order
gyrA_parC <- dplyr::filter(gyrA_parC, !grepl("ariba", ariba_ref_name))
gyrA_parC <- gyrA_parC %>% select(name, ref_ctg_change, everything())

#trim gene names in ariba_ref_name column
gyrA_parC$ariba_ref_name <- gsub(pattern = "_CP025268_K12", replacement = "", gyrA_parC$ariba_ref_name)

#filter out only mutations from our defined list
gyrA_parC <- filter(gyrA_parC, ref_ctg_change %in% gyrA_FQR | ref_ctg_change %in% parC_FQR)

#add a one or a zero in a column for gyrA and parC if a row corresponds to a FQR associated mutation
gyrA <- gyrA_parC %>% mutate(gyrA = ifelse(ariba_ref_name == "gyrA" & ref_ctg_change %in% gyrA_FQR, 1, 0)) %>% filter(ariba_ref_name == 'gyrA')
parC <- gyrA_parC %>% mutate(parC = ifelse(ariba_ref_name == "parC" & ref_ctg_change %in% parC_FQR, 1, 0)) %>% filter(ariba_ref_name == 'parC')

#produce a dataframe with a list of samples and the count of the gyrA mutations detected
gyrA_FQR_sum <- gyrA %>% select(name, gyrA)
gyrA_FQR_sum <- aggregate(gyrA_FQR_sum$gyrA, by=list(name=gyrA_FQR_sum$name), FUN=sum)

#fix the column names
colnames(gyrA_FQR_sum) <- c('name','r_gyrA*')

#produce a column that lists the specific mutations identified in gyrA
gyrA_mutations <- acast(gyrA, c(name, gyrA) ~ "gyrA_mutations", value.var= "ref_ctg_change",
      fun.aggregate=function(x) paste(x, collapse = ", ")) %>% as.data.frame()

#move sample name from rowname to column 'name' and reset the rownames
rownames(gyrA_mutations) -> gyrA_mutations$name
rownames(gyrA_mutations) <- NULL

#change column order and remove names equal to 0 or 1 - these appear due to a quirk in acast
gyrA_mutations <- gyrA_mutations %>% select(name, gyrA_mutations) %>% filter(!name %in% c(0,1))

#produce a dataframe with a list of samples and the count of the parC mutations detected
parC_FQR_sum <- parC %>% select(name, parC)
parC_FQR_sum <- aggregate(parC_FQR_sum$parC, by=list(name=parC_FQR_sum$name), FUN=sum)

#fix the column names
colnames(parC_FQR_sum) <- c('name','r_parC*')

#produce a column that lists the specific mutations identified in parC
parC_mutations <- acast(parC, c(name, parC) ~ "parC_mutations", value.var= "ref_ctg_change",
                            fun.aggregate=function(x) paste(x, collapse = ", ")) %>% as.data.frame()

#move sample name from rowname to column 'name' and reset the rownames
rownames(parC_mutations) -> parC_mutations$name
rownames(parC_mutations) <- NULL

#change column order and remove names equal to 0 or 1 - these appear due to a quirk in acast
parC_mutations <- parC_mutations %>% select(name, parC_mutations) %>% filter(!name %in% c(0,1))

#join tables
final_FQR_table <- name_list %>%
  left_join(gyrA_FQR_sum) %>% 
  left_join(parC_FQR_sum) %>%
  left_join(gyrA_mutations) %>%
  left_join(parC_mutations)

#replace NAs in columns with NA
final_FQR_table[is.na(final_FQR_table)] <- 0

#write table to csv
write_csv(final_FQR_table, "final_FQR_table.csv")
