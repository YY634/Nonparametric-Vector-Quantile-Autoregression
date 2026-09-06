######### This code converts a folder of EEG .npy files into separate R .RData files, one per subject, while extracting the subject ID from each filename.

library(RcppCNPy)

in_dir  <- "C:/Users/Desktop/NVQAR/closed_eyes"
out_dir <- "C:/Users/Desktop/NVQAR/closed_eyes_R"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
files <- sort(list.files(in_dir, pattern = "\\.npy$", full.names = TRUE))

for (f in files) {
  fname <- basename(f)
  sub_tag <- regmatches(fname, regexpr("sub-[0-9]+", fname))
  subject_id <- sub("^sub-", "", sub_tag)
  eeg_mat <- RcppCNPy::npyLoad(f)
  # check
  print(fname)
  print(dim(eeg_mat))
  print(class(eeg_mat))
  # save
  source_file <- fname
  orientation <- "channels_by_samples"
  save(eeg_mat, subject_id, source_file, orientation,
       file = file.path(out_dir, paste0("subject_", subject_id, ".RData")))
}


### closed-eye dataset is sampled at frequency fs=500 Hz, delta t = 1/500 = 0.002
### after "diluting" by use  seq(1, N, by=2), the frequency is actually 250 Hz
for (i in 1:88){
  filename <- paste0(sprintf("~/Documents/papers_until2026/Nonparametric_vector_QAR/real_data/Raw_data/closed_eye/subject_%03d", i), ".RData")
  load(filename)
  N <- dim(eeg_mat)[2]
  indseq <- seq(1, N, by=2)
  eeg <- eeg_mat[, indseq]
  save(eeg, file = filename)
  rm(list=ls())
}




############### extract the subset of data that will be input to our algorithm
sub_info <- read.table("~/Documents/papers_until2026/Nonparametric_vector_QAR/Raw_data/subjects_infor.txt", header = TRUE)
# Alzheimer's 01-36;  Control 37:65;  FTD 66-88
channels <- read.table("~/Documents/papers_until2026/Nonparametric_vector_QAR/Raw_data/eeg_channels.txt", header = FALSE)
channels <- channels[, -1]
# (F1,F2) = (row_1, row_2);   (O1,O2) = (row_9, row_10)


#### (O1, O2) in AD patients have less synchronization, less alpha rhythm, more slower oscillations, smaller amplitude, more irregular.
#### (F1, F2) in FTD patients show less alpha rhythm, more slower rhythm, less complexity, and less synchronization. Amplitude does not reduce much. 
v <- 250
 
AD <- list()
for (i in 1:36){
  filename <- paste0(sprintf("~/Documents/papers_until2026/Nonparametric_vector_QAR/Raw_data/closed_eye/subject_%03d", i), ".RData")
  load(filename)
  L <- min(100000+v, dim(eeg)[2])
  AD_tmp <- eeg[9:10, (v+1):L]
  ADT <- t(10^4 * AD_tmp) 
  colnames(ADT) <- c("O1", "O2") 
  AD[[i]] <- ADT
  rm(AD_tmp, ADT, eeg, filename)
}
rm(i)


FTD <- list()
for (i in 1:23){
  j <- i+65
  filename <- paste0("~/Documents/papers_until2026/Nonparametric_vector_QAR/Raw_data/closed_eye/subject_0", j, ".RData")
  load(filename)
  L <- min(100000+v, dim(eeg)[2])
  FTD_tmp <- eeg[1:2, (v+1):L]
  FTDD <- t(10^4*FTD_tmp)
  colnames(FTDD) <- c("F1", "F2")  
  FTD[[i]] <- FTDD
  rm(FTD_tmp, FTDD, eeg, filename)
}
rm(i)


CN_AD <- list() 
for (i in 1:29){
  j <- i+36
  filename <- paste0("~/Documents/papers_until2026/Nonparametric_vector_QAR/Raw_data/closed_eye/subject_0", j, ".RData")
  load(filename)
  L <- min(80000+v, dim(eeg)[2])
  CNT_tmp <- eeg[9:10, (v+1):L]
  CNT <- t(10^4 * CNT_tmp)
  colnames(CNT) <- c("O1", "O2") 
  CN_AD[[i]] <- CNT  
  rm(CNT, CNT_tmp, eeg, filename)
}
rm(i)


CN_FTD <- list()
for (i in 1:29){
  j <- i+36
  filename <- paste0("~/Documents/papers_until2026/Nonparametric_vector_QAR/Raw_data/closed_eye/subject_0", j, ".RData")
  load(filename)
  L <- min(80000+v, dim(eeg)[2])
  CNT_tmp <- eeg[1:2, (v+1):L]
  CNT <- t(10^4 * CNT_tmp)
  colnames(CNT) <- c("F1", "F2") 
  CN_FTD[[i]] <- CNT
  rm(CNT, CNT_tmp, eeg, filename)
}

save(AD, FTD, CN_AD, CN_FTD, file = "~/Documents/papers_until2026/Nonparametric_vector_QAR/JASA_revision_round2/real_data/input.RData")




