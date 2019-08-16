library(readr)
df <- read_table2("SNP_matrix_ST744.txt", skip = 1)

compare <- function(x){
  i <- 4
  for(i in 4:(ncol(x)-4)){
    print(sum(x[,3] != x[,i]))
    }
  }

compare(df)