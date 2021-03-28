setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Combine all annotations
files = list.files(pattern = "aaAnnotation.csv")
df = NULL
for (i in 1:length(files)) {
  if (i == 1) {
    df = read.csv(files[i])
  }
  else {
    temp = read.csv(files[i])
    df = merge(df, temp, by = "sequence")
  }
}

# Label irrelavent mutations as "X"
mut = c("A222V", "N439K", "S477N", "D614G")
for (i in 1:length(mut)) {
  col = substr(mut[i], 1, 4)
  wt = substr(mut[i], 1, 1)
  mu = substr(mut[i], 5, 5)
  df[!(df[,col] %in% c(wt, mu)),i+1] = "X"
}

write.csv(df, file = "combined_aaAnnotation.csv", row.names = FALSE)

ID = c("EPI_ISL_462861", 
       "hCoV-19/Australia/VIC122/2020|EPI_ISL_419724|2020-03-18", 
       "hCoV-19/Australia/VIC2215/2020|EPI_ISL_521955|2020-06-14")
df_test = df[df$sequence %in% ID,]
write.csv(df_test, file = "testing_aaAnnotation.csv", row.names = FALSE)
