library("R.utils")

args = commandArgs(trailingOnly=TRUE)

files <- list.files(path = ".", pattern = glob2rx(paste0(args[1],"_*ldist.csv")))

list_lines=c()
DF <-  read.csv(files[1],header=F)
list_lines = c(list_lines,countLines(files[1]))


for (f in files[-1]){
  df <- read.csv(f,header=F)      # read the file
  DF <- rbind(DF, df)    # append the current file
  list_lines = c (list_lines,countLines(f))
}
#writing the appended file  
write.csv(DF, paste0(args[1],".all_ldist_concat.csv"), 
          row.names=FALSE, quote=FALSE)

freq=as.data.frame(table(DF[2]))

colnames(freq) <- c("SNP","count")

tot = sum(freq$count)
freq$Frequency <- freq$count/tot


write.csv(freq, paste0(args[1],".tbl_freq.csv"), 
          row.names=FALSE, quote=FALSE)
               
               #
               
               
tab <- matrix(c("nsim", length(files), "mean_samples",mean(list_lines),
                "median_samples", median(list_lines), "max_samples", 
                max(list_lines), "min_samples", min(list_lines)), 
              ncol=2, byrow=TRUE)
tab <- as.table(tab)

write.table(tab, file = paste0(args[1],".info"),quote = FALSE,
            row.names = FALSE, col.names = FALSE)

tab <- cbind(list_lines)

write.table(tab, file = paste0(args[1],".N_count"),quote = FALSE,
            row.names = FALSE, col.names = FALSE)




files <- list.files(path = ".", pattern = glob2rx(paste0(args[1],"_*.cl_r")))

DF <-  read.table(files[1],header=F)

for (f in files[-1]){
  df <- read.table(f,header=F)      # read the file
  DF <- rbind(DF, df)    # append the current file
}
#writing the appended file  
write.table(DF, paste0(args[1],".all_cl_r_concat"), 
          row.names=FALSE, quote=FALSE,col.names=FALSE)

#files <- list.files(path = ".", pattern = glob2rx(paste0(args[1],"_*.clnn_r")))

#DF <-  read.table(files[1],header=F)

#for (f in files[-1]){
#  df <- read.table(f,header=F)      # read the file
#  DF <- rbind(DF, df)    # append the current file
#}
#writing the appended file  
#write.table(DF, paste0(args[1],".all_clnn_r_concat"), 
#          row.names=FALSE, quote=FALSE,col.names=FALSE)
