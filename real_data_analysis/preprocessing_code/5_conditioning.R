


load("~/Documents/papers_until2026/Nonparametric_vector_QAR/JASA_revision_round2/real_data/represents.RData")

tims <- c(15697, 20643, 25036, 32316, 41007, 50477, 58984, 60520, 69097, 76045, 81452, 93170) 
conds_AD <- AD_medoid[tims,]; rm(tims)

tims <- c(4075, 10199, 23990, 33443, 40025, 53863, 60100, 61641,  63000,  69124, 71917, 77081)  
conds_CNAD <- CNAD_medoid[tims,];  rm(tims)

tims <- c(1067, 19998, 25193, 33525, 44845, 44848,  52982, 60475, 60546, 61013, 70002, 80145)  
conds_FTD <- FTD_medoid[tims,];   rm(tims)

tims <- c(4163, 7918, 14956, 18220, 20352, 44920, 50031, 54001, 56430, 66616, 69527, 76237)
conds_CNFTD <- CNFTD_medoid[tims,];  rm(tims)


save(conds_AD, conds_CNAD, conds_FTD, conds_CNFTD, file = "~/Documents/papers_until2026/Nonparametric_vector_QAR/JASA_revision_round2/real_data/conds.RData")


