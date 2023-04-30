library(sas7bdat)
mydata <- read.sas7bdat("merge_data2_1.sas7bdat", encoding="", debug=FALSE)

# VE
dataset_1 <- mydata[mydata$arm == "Vitamin E" | mydata$arm == "Placebo" ,]
# Donepezil
dataset_2 <- mydata[mydata$arm == "Donepezil" | mydata$arm == "Placebo",]


# VE
id_1 <- dataset_1$Patient_ID[dataset_1$viscode == "bl"]
# Donepezil
id_2 <- dataset_2$Patient_ID[dataset_2$viscode == "bl"]

# VE
dataset_1 <- dataset_1[dataset_1$Patient_ID %in% id_1,]
# Donepezil
dataset_2 <- dataset_2[dataset_2$Patient_ID %in% id_2,]

getmode <- function(x)
{
  if (length(x) > 0){
    return(as.numeric(names(table(x))[table(x) == max(table(x))]))
  }
  else{
    return(0)
  }
}

# Donepezil
theshold <- 18
baseline <- dataset_2$cdsum[dataset_2$viscode == "bl"]
y1 <- c()
y2 <- c()
y3 <- c()
y4 <- c()
y5 <- c()
y6 <- c()
gender <- c()
age <- c()
genotype <- c()
cototscr <- c()
adtotscr <- c()
evgds <- c()
mmscore <- c()
trt <- c()
for (i in 1:length(id_2)) {
  id_cur <- id_2[i]
  y1_cur <- dataset_2$cdsum[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "m06"]
  if (length(y1_cur) == 0){
    y1[i] <- NA
  }
  else if (y1_cur < 0 | y1_cur > theshold){
    y1[i] <- NA
  }
  else{
    y1[i] <- y1_cur
  }
  y2_cur <- dataset_2$cdsum[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "m12"]
  if (length(y2_cur) == 0){
    y2[i] <- NA
  }
  else if (y2_cur < 0 | y2_cur > theshold){
    y2[i] <- NA
  }
  else{
    y2[i] <- y2_cur
  }
  y3_cur <- dataset_2$cdsum[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "m18"]
  if (length(y3_cur) == 0){
    y3[i] <- NA
  }
  else if (y3_cur < 0 | y3_cur > theshold){
    y3[i] <- NA
  }
  else{
    y3[i] <- y3_cur
  }
  y4_cur <- dataset_2$cdsum[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "m24"]
  if (length(y4_cur) == 0){
    y4[i] <- NA
  }
  else if (y4_cur < 0 | y4_cur > theshold){
    y4[i] <- NA
  }
  else{
    y4[i] <- y4_cur
  }
  y5_cur <- dataset_2$cdsum[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "m30"]
  if (length(y5_cur) == 0){
    y5[i] <- NA
  }
  else if (y5_cur < 0 | y5_cur > theshold){
    y5[i] <- NA
  }
  else{
    y5[i] <- y5_cur
  }
  y6_cur <- dataset_2$cdsum[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "m36"]
  if (length(y6_cur) == 0){
    y6[i] <- NA
  }
  else if (y6_cur < 0 | y6_cur > theshold){
    y6[i] <- NA
  }
  else{
    y6[i] <- y6_cur
  }
  gender[i] <- dataset_2$ptgender[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  age[i] <- dataset_2$infage6[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  genotype[i] <- dataset_2$genotype2[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  cototscr[i] <- dataset_2$cototscr[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  adtotscr[i] <- dataset_2$adtotscr[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  evgds[i] <- dataset_2$evgds[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  mmscore[i] <- dataset_2$mmscore[dataset_2$Patient_ID == id_cur & dataset_2$viscode == "bl"]
  trt[i] <- dataset_2$arm[dataset_2$Patient_ID == id_cur & dataset_2$viscode == 'bl'] != "Placebo"
}
e1 <- y1 - baseline
e2 <- y2 - baseline
e3 <- y3 - baseline
e4 <- y4 - baseline
e5 <- y5 - baseline
e6 <- y6 - baseline

age_del <- age < 0
gender_del <- gender < 0
mmscore_del <- mmscore < 0
cototscr_del <- cototscr < 0
adtotscr_del <- adtotscr < 0
evgds_del <- evgds < 0
del_term <- age_del | gender_del | mmscore_del | cototscr_del | adtotscr_del | evgds_del

e1 <- e1[!del_term]
e2 <- e2[!del_term]
e3 <- e3[!del_term]
e4 <- e4[!del_term]
e5 <- e5[!del_term]
e6 <- e6[!del_term]
age <- age[!del_term]
baseline <- baseline[!del_term]
gender <- gender[!del_term]
mmscore <- mmscore[!del_term]
cototscr <- cototscr[!del_term]
adtotscr <- adtotscr[!del_term]
evgds <- evgds[!del_term]
genotype <- genotype[!del_term]
genotype <- genotype - 3
trt <- trt[!del_term]
trt <- as.numeric(trt)
N <- length(genotype)


wideData_2 <- data.frame(id = 1:N, trt = trt, e1 = e1, e2 = e2, e3 = e3, e4 = e4,
                         e5 = e5, e6 = e6, BL = baseline, AGE = age, GENDER = as.factor(gender),
                         MMSCORE = mmscore,COTOTSCR = cototscr, ADTOTSCR = adtotscr,
                         EVGDS = evgds, GENOTYPE = as.factor(genotype))

cor_dataset <- wideData_2[,c(9,12,13,14,15)]
cor(cor_dataset)

placebo_group <- wideData_2$e6[wideData_2$trt == 0]
treatment_group <- wideData_2$e6[wideData_2$trt == 1]
t.test(placebo_group,treatment_group)


t1 <- Sys.time()
set.seed(9)
tr <- MMRM_tree(wideData_2,e~factor(Z)*factor(trt)*factor(time)+us(time|id),
                c("BL","AGE","GENDER","MMSCORE","COTOTSCR","ADTOTSCR","EVGDS","GENOTYPE"),
                Prunemethod = "Bootstrap",level_conti = 30,min_N = 30,p.value = 0.1,lambda = c(0, 2, 3, 4, log(378)))
t2 <- Sys.time()
t2-t1

print(tr)

subgroup_1 <- wideData_2[wideData_2$MMSCORE <= 24,]
dim(subgroup_1)
sum(subgroup_1$trt == 1)
sum(subgroup_1$trt == 0)
subgroup_2 <- wideData_2[wideData_2$MMSCORE > 24 & wideData_2$ADTOTSCR <= 43,]
dim(subgroup_2)
sum(subgroup_2$trt == 1)
sum(subgroup_2$trt == 0)
subgroup_3 <- wideData_2[wideData_2$MMSCORE > 24 & wideData_2$ADTOTSCR > 43 & wideData_2$COTOTSCR <= 15,]
dim(subgroup_3)
sum(subgroup_3$trt == 1)
sum(subgroup_3$trt == 0)
subgroup_4 <- wideData_2[wideData_2$MMSCORE > 24 & wideData_2$ADTOTSCR > 43 & wideData_2$COTOTSCR > 15,]
dim(subgroup_4)
sum(subgroup_4$trt == 1)
sum(subgroup_4$trt == 0)
subgroup_combine <- rbind(subgroup_1,subgroup_2)

node <- subgroup_1
#node <- subgroup_2
#node <- subgroup_3
#node <- subgroup_4
#node <- subgroup_combine

longData_cur <- reshape(node, varying=grep(pattern = "e",names(node),value = TRUE),
                        direction="long", sep="", idvar="id")
longData_cur <- longData_cur[order(longData_cur$trt,longData_cur$id,longData_cur$time),]
longData_cur$time <- factor(longData_cur$time)
longData_cur$id <- factor(longData_cur$id)
mmrm_cur <- mmrm(e~factor(trt)*factor(time)+us(time|id),data = longData_cur)
contrast <- numeric(length(mmrm_cur$beta_est))
contrast[2] <- contrast[12] <- 1
df_1d(mmrm_cur, contrast)
