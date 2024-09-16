

data(aaindex, package = "seqinr")

idx <- sample(1:length(aaindex), 5)

for (i in idx) {
  print(aaindex[[i]]$D)
}

prompt <- "
I am evaluating amino acid indices from a dataset that lists various amino acid properties. The evaluation should determine whether the index is suitable to quantify pairwise amino acid  \"similarity\" if the similarity is to be used to evaluate whether the structure of the genetic code makes it robust against the effect of point mutations. Therefore it is important to exclude indices which would lead to circular reasoning, i.e. where the reported values derive from the structure of the code, such as codon adjacency, or the unequal frequency of amino acids in the code. I need help categorizing a specific index based on its origin and potential biases. Here are five descriptions indices from the table:

\"[Insert Index Description Here]\"

Based on this description, would you categorize this index as:
(A). Lab-derived: Properties determined from isolated amino acids in controlled experiments.
(B). Behavioral Observations: Properties based on observed structural propensities of amino acids in proteins, or other roles in metabolism, that might contain selection bias based on codon adjancency or different frequencies of representation of amino acids in the code.
(C). Code-derived: Properties that entail from the features or biases inherent in the genetic code itself.

Please provide a category recommendation along with a concise reasoning. Additionally, please estimate your confidence level (as a percentage) in the categorization decision? If the confidence is low, please specify any additional information or context that might help refine the assessment. You may recommend none, one, two, or all three of the categories with their respective confidence levels.
"


# (Evolution of the code -> symmetry breaking)
#
cat(sprintf("\n%s",sample(c(AADAT$nam[-21], sample(AADAT$nam[-21],8)))))

x <- t(as.matrix(aaindex[[1]]$I))

# produce a dataframe
aaAll <- as.data.frame(t(aaindex[[1]]$I))
colnames(aaAll) <- names(aaindex[[1]]$I)
rownames(aaAll)[1] <-  aaindex[[1]]$H

iNA <- logical(length(aaindex))
iNA[1] <- any(is.na(aaindex[[1]]$I))

for (i in 2:length(aaindex)) {
  aaAll <- rbind(aaAll, aaindex[[i]]$I)
  rownames(aaAll)[i] <-  aaindex[[i]]$H
  iNA[i] <- any(is.na(aaindex[[i]]$I))
}

which(iNA)

for (i in 1:length(which(iNA))) {
  idx <- which(iNA)[i]
  print(aaindex[[idx]]$H)
  print(aaindex[[idx]]$D)
  print(aaindex[[idx]]$T)
  print("===============================")
}

aain


# standardize

aaStd <- as.data.frame(t(scale(t(aaAll), center=TRUE, scale=TRUE)))
head(aaStd)

row_means <- apply(aaStd, 1, mean)
row_sds <- apply(aaStd, 1, sd)
summary(row_means)
summary(row_sds)

aaStd <- aaStd[! is.na(row_means), ]

pcaAA <- princomp(aaStd)
pcaAA

plot(pcaAA)

highest_correlation <- function(pca = pcaAA, idx = 1) {
  # Get the specified component's loadings
  load <- pca$loadings[, idx]

  # Find the index (row) with the highest absolute correlation
  index <- which.max(abs(loadings_vec))

  return(index)
}

iInclude <- ! iNA

Curators <- strsplit("
Francesca
Lara
Rachel
Raihan
Ariuntuya
Nabiha
James
Jacob
Claire
Ria
Kai
Minsoo
Jaycee
Batool
Simon
Ian
Kinley
Anish
Lauren
Zehao
Quinn
Ayesha
Suhani
Evan
Tohya
Harry
Peter
Dechen
Boris", "\\s+")[[1]]

Curators <- Curators[Curators != ""]

allNames <- rep(Curators, 19)

nRow <- length(aaindex) - length(which(iNA))
curationTable <- data.frame(IDX        = integer(nRow),
                           ID          = character(nRow),
                           Curator     =  character(nRow),
                           AI.prop     = numeric(nRow),
                           AI.behave   = numeric(nRow),
                           AI.derived  = numeric(nRow),
                           Cur.prop    = numeric(nRow),
                           Cur.behave  = numeric(nRow),
                           Cur.derived = numeric(nRow),
                           Description = character(nRow),
                           Title       = character(nRow),
                           Notes       = character(nRow))
for (i in 1:20) {
  curationTable <- cbind(curationTable, numeric(nRow))
}
colnames(curationTable)[13:32] <- c("Ala", "Cys", "Asp", "Glu", "Phe",
                                    "Gly", "His", "Ile", "Lys", "Leu",
                                    "Met", "Asn", "Pro", "Gln", "Arg",
                                    "Ser", "Thr", "Val", "Trp", "Tyr")
iRows <- (1:length(aaindex))[! iNA]
tail(iRows)

for (i in 1:nRow) {
  iAA <- iRows[i]
  curationTable[i, "IDX"]         <- iAA
  curationTable[i, "ID"]          <- aaindex[[iAA]]$H
  curationTable[i, "Curator"]     <- allNames[i]
  curationTable[i, "AI.prop"]     <- 0
  curationTable[i, "AI.behave"]   <- 0
  curationTable[i, "AI.derived"]  <- 0
  curationTable[i, "Cur.prop"]    <- 0
  curationTable[i, "Cur.behave"]  <- 0
  curationTable[i, "Cur.derived"] <- 0
  curationTable[i, "Description"] <- aaindex[[iAA]]$D
  curationTable[i, "Title"]       <- aaindex[[iAA]]$T
  curationTable[i, "Notes"]       <- ""
  for (Aaa in colnames(curationTable)[13:32]) {
    curationTable[i, Aaa] <- aaindex[[iAA]]$I[Aaa]
  }
}

# sample() now - everyone is already present the most number of times
curationTable$Curator <- sample(curationTable$Curator)

write.csv(curationTable, file = "curationTable.csv", row.names = FALSE)

# ===========

sds <- pcaAA$sdev
vari <- sds^2
variTot <- sum(vari)
variCum <- cumsum(vari)
variProp <- vari/variTot
variPropCum <- variCum / variTot
plot(variPropCum, type="l", col="#AA0000")
abline(h=0.9, col="#AABBFF")
barplot(variPropCum)

# Get the first PC

myLoad <- pcaAA$loadings
myLoad[,1]
plot(myLoad[,1], myLoad[,2], type = "n")
text(myLoad[,1], myLoad[,2], labels = rownames(myLoad))

plot(myLoad[,2], myLoad[,3], type = "n")
text(myLoad[,2], myLoad[,3], labels = rownames(myLoad))


cor1 <- numeric()
cor2 <- numeric()
cor3 <- numeric()

myLoad[,1]
aaindex[[77]]$I
cor(myLoad[,1], aaindex[[77]]$I)

for (i in 1:length(aaindex)) {
  cor1[i] <- cor(myLoad[,1], aaindex[[i]]$I)
  cor2[i] <- cor(myLoad[,2], aaindex[[i]]$I)
  cor3[i] <- cor(myLoad[,3], aaindex[[i]]$I)
}

summary(cor1)
summary(cor2)
summary(cor3)

aaindex[[which(cor1 == max(cor1, na.rm = TRUE))]] # 211
aaindex[[which(cor2 == max(cor2, na.rm = TRUE))]] # 81
aaindex[[which(cor3 == max(cor3, na.rm = TRUE))]] # 470

# === Reading from a Google sheet =========
tmp <- tmp[!is.na(tmp$Cur.derived), ]
tmp <- tmp[tmp$Cur.derived == 0, ]

colnames(tmp)[14:33]


aaStd <- as.data.frame(t(scale(t(tmp[14:33]), center=TRUE, scale=TRUE)))
head(aaStd)

row_means <- apply(aaStd, 1, mean)
row_sds <- apply(aaStd, 1, sd)
summary(row_means)
summary(row_sds)
rownames(aaStd) <- tmp$ID

aaStd <- aaStd[! is.na(row_means), ]

pcaAA <- princomp(aaStd)
print(pcaAA)


sds <- pcaAA$sdev
vari <- sds^2
variTot <- sum(vari)
variCum <- cumsum(vari)
variProp <- vari/variTot
variPropCum <- variCum / variTot
barplot(variPropCum)
abline(h=0.9, col="#AABBFF")

# Get the first PC

myLoad <- pcaAA$loadings
summary(myLoad[,1])
summary(myLoad[,2])
summary(myLoad[,3])


myLoad[,1]
plot(myLoad[,1], myLoad[,2], type = "n")
text(myLoad[,1], myLoad[,2], labels = rownames(myLoad))

plot(myLoad[,2], myLoad[,3], type = "n")
text(myLoad[,2], myLoad[,3], labels = rownames(myLoad))


cor1 <- numeric()
cor2 <- numeric()
cor3 <- numeric()


cor(myLoad[,1], t(aaStd[77, ]))

for (i in 1:nrow(aaStd)) {
  cor1[i] <- cor(myLoad[,1], t(aaStd[i, ]))
  cor2[i] <- cor(myLoad[,2], t(aaStd[i, ]))
  cor3[i] <- cor(myLoad[,3], t(aaStd[i, ]))
}

summary(cor1)
summary(cor2)
summary(cor3)

rownames(aaStd)[which(cor1 == max(cor1, na.rm = TRUE))] # "NISK860101"

rownames(aaStd)[which(cor2 == max(cor2, na.rm = TRUE))] # "RADA880105"
rownames(aaStd)[which(cor2 == min(cor2, na.rm = TRUE))] # "FAUJ880104"

rownames(aaStd)[which(cor3 == max(cor3, na.rm = TRUE))] # "KUMS000104"

plot(myLoad[,1], aaStd["NISK860101", ])
plot(myLoad[,2], aaStd["FAUJ880104", ])

# END
