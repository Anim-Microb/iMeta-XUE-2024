#################################PCA####################
library(dplyr)
library(ggplot2)
library(vegan)
data = read.csv(file.choose(), row.names = 1)
  as.data.frame() %>%
  mutate(group = rep(c('LFCR','HFCR')))
pca = prcomp(data[,1:ncol(data)-1],scale. = TRUE)
pca.var = pca$sdev^2 
pca.var.per = round(pca.var/sum(pca.var)*100,2)
data.frame(PC = as.character(paste('PC',1:nrow(data),sep = '')),
           value = pca.var.per) %>%
  ggplot(aes(PC,value))+
  geom_bar(stat = 'identity', fill = 'white', color = 'black')+
  geom_hline(yintercept = pca.var.per[1]*1.1, color = 'white')+
  scale_x_discrete(limits = c(paste('PC',1:nrow(data),sep = '')))+
  theme_classic()+
  scale_y_continuous(expand = c(0,0))+
  geom_text(aes(y = value + 1,label = paste(value,'%',sep = '')),size = 2.5)+
  labs(x = '',y = '',title = 'ScreePlot')+
  theme(axis.text = element_text(size = 11,color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5))
as.data.frame(pca$x) %>%
  mutate(group = data$group) %>%
  ggplot(aes(PC1,PC2,color = group))+
  geom_point(size = 2)+
  theme_classic()+
  labs(x = paste('PC1(',pca.var.per[1],'%)',sep = ''),
       y = paste('PC2(',pca.var.per[2],'%)',sep = ''))+
  stat_ellipse(level = 0.68)+
  theme(axis.text = element_text(size = 11,color = 'black'),
        axis.ticks = element_line(color = 'black'),
        plot.title = element_text(hjust = 0.5))

#####################volcano plot#####################
library(ggplot2)
log <- read.csv(file.choose(), row.names = 1)
log2Foldchange <- read.csv(file.choose(), row.names = 1)
log10pvalue <- read.csv(file.choose(), row.names = 1)
expression <- ifelse(log10pvalue >= 1.5 & log2Foldchange < 0, "Down-regulated", ifelse(log10pvalue >= 1.5 & log2Foldchange > 0, "up-regulated"))
ggplot(log, aes(log2Foldchange,log10pvalue)) + geom_point(aes(ccolor = expression)) + scale_y_continuous(limits = c(0,1.5))
#################PERMANOVA-MTT-MTB-MTP#########################
mtb <- read.csv(file.choose(), row.names = 1)
library(PERMANOVA)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(vegan)
d <- dist(mtb, method = "euclidean")
mtb$group <- as.factor(mtb$group)
adonis2(d ~ group, data = mtb, permutations = 999, method = 'bray')

mtt <- read.csv(file.choose(),row.names = 1)
d <- dist(mtt, method = "euclidean")
mtt$group <- as.factor(mtt$group)
adonis2(d ~ group, data = mtt, permutations = 999, method = 'bray')

mtp <- read.csv(file.choose(),row.names = 1)
d <- dist(mtp, method = "euclidean")
mtp$group <- as.factor(mtp$group)
adonis2(d ~ group, data = mtp, permutations = 999, method = 'bray')

#######################TwoSampleMR#########################
install.packages("TwoSampleMR-0.4.22.tar")
system.file
install.packages("gsmr_1.0.8.tar")
install.packages("MendelianRandomization_0.7.0.tar")
install.packages(c('survey'));
remotes::install_github("MRCIEU/TwoSampleMR")
install.packages("TwoSampleMR-0.5.6.tar", repos = NULL, type = "source")
install.packages("TwoSampleMR-0.5.6.tar")
install.packages("ieugwasr")-
  install.packages("markdown")
install.packages("MendelianRandomization")
install.packages("mr.raps")
install.packages("randomForest")
install.packages("meta")
install.packages("MRPRESSO")-
  install.packages("MRInstruments")-
  install.packages("RadialMR")-
  install.packages("MRMix")-
  install.packages("glmnet")
install.packages("pbapply")
library(TwoSampleMR)
exp_dat <- read_exposure_data(
  file.choose(),
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "a1_freq",
  pval_col = "P",
  phenotype_col = "id.exposure"
)
out_dat1 <- read_outcome_data(
  file.choose(),
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.outcome",
  se_col = "se.outcome",
  effect_allele_col = "effect_allele.outcome",
  other_allele_col = "other_allele.outcome",
  eaf_col = "a1_freq",
  pval_col = "pval.outcome",
  phenotype_col = "id.outcome"
)
dat1 <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat1
)
res1 <- mr(dat1)
out_dat2 <- read_outcome_data(file.choose(),sep = ",",snp_col = "SNP",beta_col = "beta.outcome",se_col = "se.outcome",effect_allele_col = "effect_allele.outcome",other_allele_col = "other_allele.outcome",eaf_col = "a1_freq",pval_col = "pval.outcome",phenotype_col = "id.outcome")
dat2 <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat2)
res2 <- mr(dat2)
cbind()
rbind()
for (i in 1:259) {out_dat[i] <- read.csv(paste0("assoc_results",i, ".qassoc.csv"), header = TRUE, sep = ";")}
temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)
paste0("assoc_results",i, ".qassoc.csv")
file_names <- paste0("out_dat_", 1:259, ".csv")
out_dat2 <- read_outcome_data(paste0("assoc_results",i, ".qassoc.csv"),sep = ",",snp_col = "SNP",beta_col = "beta.outcome",se_col = "se.outcome",effect_allele_col = "effect_allele.outcome",other_allele_col = "other_allele.outcome",eaf_col = "a1_freq",pval_col = "pval.outcome",phenotype_col = "id.outcome")
file_names <- paste0("out_dat_", 1:259, ".csv")
out_dat_list <- lapply(file_names, function(file_name) {
  file <- paste0(file_name)
  data <- read.csv(file)
  data <- clean_data(data)
  return(data$out_dat)
})
library(TwoSampleMR)
file_path <- "/Users/xmyzyf/Desktop/onesampleMR/outcome/"
file_names <- paste0("assoc_results", 1:259, ".qassoc.csv")
out_dat_list <- list()
for (i in 1:length(file_names)) {
  file <- paste0(file_path, file_names[i])
  data <- read_outcome_data(file, sep = ",",
                            snp_col = "SNP",
                            beta_col = "beta.outcome",
                            se_col = "se.outcome",
                            effect_allele_col = "effect_allele.outcome",
                            other_allele_col = "other_allele.outcome",
                            eaf_col = "a1_freq",
                            pval_col = "pval.outcome",
                            phenotype_col = "id.outcome")
  out_dat_list[[i]] <- data
}
out_dat <- do.call(rbind, out_dat_list)
dat_list <- lapply(out_dat_list, function(out_dat) {
  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  res <- mr(dat)
  return(res)
})


write.csv(dat_list, file = "res1")
write.csv(dat_list, file = "res1-1")
fit <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))


######################fit#########################
data5 <- read.csv(file.choose(), row.names = 1)
ggplot(data = data5,
       mapping = aes( x = parity, y= ace, color = breed)) + geom_smooth() #两条线
ggplot(data = cagcor,
       mapping = aes( x = PI, y= CAG4)) + geom_point() + stat_smooth(method = "lm") #直线
ggplot(data = data5,
       mapping = aes( x = parity, y= ace, color = parity)) + geom_point(aes(color = parity)) + stat_smooth() 
#
ggplot(data = mtcars, mapping = aes(x = wt, y = mpg)) + geom_point(aes(color = am)) + stat_smooth()
#https://blog.csdn.net/ethmery/article/details/112583698
pyruvate <- read.csv(file.choose(), row.names = 1)
propionate <- read.csv(file.choose(), row.names = 1)
mcp <- read.csv(file.choose(), row.names = 1)
chea <- read.csv(file.choose(), row.names = 1)
chea$Selenomonas.bovis 
chea$Pyruvate
library(ggplot2)
library(Hmisc)
lmpyruvate <- ggplot(data = pyruvate,
                     mapping = aes( x = Selenomonas.bovis, y= Pyruvate)) + geom_point() + stat_smooth(method = "lm") 
lmpyruvate
pyruvate <- as.matrix(pyruvate)
corpyruvate <- rcorr(pyruvate, type = "spearman")
corpyruvate$r
corpyruvate$P
lmpropionate <- ggplot(data = propionate,
                       mapping = aes( x = Selenomonas.bovis, y= Propionate)) + geom_point() + stat_smooth(method = "lm") 
lmpropionate
propionate <- as.matrix(propionate)
corpropionate <- rcorr(propionate, type = "spearman")
corpropionate$r
corpropionate$P
lmmcp <- ggplot(data = mcp,
                mapping = aes( x = Selenomonas.bovis, y= mcp)) + geom_point() + stat_smooth(method = "lm") 
lmmcp
mcp <- as.matrix(mcp)
cormcp <- rcorr(mcp, type = "spearman")
cormcp$r
cormcp$P
lmchea <- ggplot(data = chea,
                 mapping = aes( x = Selenomonas.bovis, y= cheA)) + geom_point() + stat_smooth(method = "lm") 
lmchea
chea <- as.matrix(chea)
corchea <- rcorr(chea, type = "spearman")
corchea$r
corchea$P
###########################median analysis###################
install.packages("mediation")
library(mediation)
data <- read.csv(file.choose(), row.names = 1)
data$efficiency
data$S.bovis
data$propionate
summary(model)
m <- lm(propionate ~ bacteria , data = data)
y <- lm(efficiency ~ bacteria + propionate , data = data)
m
y
summary(m)
summary(y)
model <- mediate(m, y, sims = 50, treat = "bacteria",  mediator = "propionate", boot = TRUE)
summary(model)
plot(model)