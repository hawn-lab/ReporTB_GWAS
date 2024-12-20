library(tidyverse)
library(broom)
library(matrixTests)

#### Data ####
attach("../metadata/ReportTB_metadata_clean.RData")

meta_FULL <- meta_FULL %>% mutate(smokhx2 = recode(smokhx, "past"="current"))
meta_CC <- meta_CC %>% mutate(smokhx2 = recode(smokhx, "past"="current"))

#### Full cohort ####
summary(meta_FULL$bl_age)
sd(meta_FULL$bl_age,na.rm = TRUE)

summary(meta_FULL$bmi)
sd(meta_FULL$bmi,na.rm = TRUE)

meta_FULL %>% count(sex) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% count(race) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% count(diabetes) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% count(smokhx) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% count(bl_hiv) %>% mutate(pct = n/sum(n)*100)

meta_FULL %>% filter(bl_hiv=="Y") %>% count(bl_art) %>% mutate(pct = n/sum(n)*100)
summary(meta_FULL$lab_cd4)
sd(meta_FULL$lab_cd4,na.rm = TRUE)

#### Full cohort by group ####
table(meta_FULL$tb)

meta_FULL %>% 
  group_by(tb) %>% 
  summarise(age = mean(bl_age), age.s = sd(bl_age),
            age.min = min(bl_age), age.max= max(bl_age))

summary(meta_FULL[meta_FULL$tb=="contact",]$bmi)
sd(meta_FULL[meta_FULL$tb=="contact",]$bmi,na.rm = TRUE)
summary(meta_FULL[meta_FULL$tb=="culture_confirm",]$bmi)
sd(meta_FULL[meta_FULL$tb=="culture_confirm",]$bmi,na.rm = TRUE)

meta_FULL %>% group_by(tb) %>% count(sex) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% group_by(tb) %>% count(race) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% group_by(tb) %>% count(smokhx2) %>% mutate(pct = n/sum(n)*100)

meta_FULL %>% group_by(tb) %>% count(diabetes) %>% mutate(pct = n/sum(n)*100)
# summary(meta_FULL[meta_FULL$tb=="contact",]$lab_hgapct)
# sd(meta_FULL[meta_FULL$tb=="contact",]$lab_hgapct,na.rm = TRUE)
summary(meta_FULL[meta_FULL$tb=="culture_confirm",]$lab_hgapct)
sd(meta_FULL[meta_FULL$tb=="culture_confirm",]$lab_hgapct,na.rm = TRUE)

meta_FULL %>% group_by(tb) %>% count(bl_hiv) %>% mutate(pct = n/sum(n)*100)
meta_FULL %>% filter(bl_hiv=="Y") %>% 
  group_by(tb) %>% count(bl_art) %>% mutate(pct = n/sum(n)*100)
summary(meta_FULL[meta_FULL$tb=="contact",]$lab_cd4)
sd(meta_FULL[meta_FULL$tb=="contact",]$lab_cd4,na.rm = TRUE)
summary(meta_FULL[meta_FULL$tb=="culture_confirm",]$lab_cd4)
sd(meta_FULL[meta_FULL$tb=="culture_confirm",]$lab_cd4,na.rm = TRUE)

# Under 5
meta_FULL %>% filter(bl_age <5) %>% count(tb)

# ltbi treat
meta_FULL %>% filter(max_ltbi_treat_days>0 & max_ltbi_treat_days<Inf) %>% group_by(tb) %>% summarise(n=n(),m=mean(max_ltbi_treat_days), sd=sd(max_ltbi_treat_days))

#### Case-control by group ####
table(meta_CC$tb)

meta_CC %>% 
  group_by(tb) %>% 
  summarise(age = mean(bl_age), age.s = sd(bl_age),
            age.min = min(bl_age), age.max= max(bl_age))

summary(meta_CC[meta_CC$tb=="contact",]$bmi)
sd(meta_CC[meta_CC$tb=="contact",]$bmi,na.rm = TRUE)
summary(meta_CC[meta_CC$tb=="culture_confirm",]$bmi)
sd(meta_CC[meta_CC$tb=="culture_confirm",]$bmi,na.rm = TRUE)

meta_CC %>% group_by(tb) %>% count(sex) %>% mutate(pct = n/sum(n)*100)
meta_CC %>% group_by(tb) %>% count(race) %>% mutate(pct = n/sum(n)*100)
meta_CC %>% group_by(tb) %>% count(smokhx2) %>% mutate(pct = n/sum(n)*100)
meta_CC %>% group_by(tb) %>% count(bl_hiv) %>% mutate(pct = n/sum(n)*100)
meta_CC %>% group_by(tb) %>% count(diabetes) %>% mutate(pct = n/sum(n)*100)

# Under 5
meta_CC %>% filter(bl_age <5) %>% count(tb)

# ltbi treat
meta_CC %>% filter(max_ltbi_treat_days>0 & max_ltbi_treat_days<Inf) %>% group_by(tb) %>% summarise(n=n(),m=mean(max_ltbi_treat_days), sd=sd(max_ltbi_treat_days))


#### sceQTL ####
attach("../ReportTB_sceQTL/data/scRNAseq/ReportTB_voom.RData")

meta_SC <- voom_B$targets %>% 
  distinct(sample_id,tb, bl_age, bmi, sex, race, smokhx2,bl_hiv, diabetes)

table(meta_SC$tb)

meta_SC %>% 
  summarise(age = mean(bl_age), age.s = sd(bl_age),
            age.min = min(bl_age), age.max= max(bl_age))

summary(meta_SC$bmi)
sd(meta_SC$bmi,na.rm = TRUE)

meta_SC %>% count(sex) %>% mutate(pct = n/sum(n)*100)
meta_SC %>% count(race) %>% mutate(pct = n/sum(n)*100)
meta_SC %>% count(smokhx2) %>% mutate(pct = n/sum(n)*100)
meta_SC %>% count(bl_hiv) %>% mutate(pct = n/sum(n)*100)
meta_SC %>% count(diabetes) %>% mutate(pct = n/sum(n)*100)

#### Statistics: full cohort ####
wilcox.x <- meta_FULL %>% 
  filter(tb == "contact") %>% 
  select(bl_age, bmi, lab_cd4) 
wilcox.y <- meta_FULL %>% 
  filter(tb == "culture_confirm") %>% 
  select(bl_age, bmi, lab_cd4) 

#Numeric
col_wilcoxon_twosample(wilcox.x, wilcox.y, alternative = "two.sided",
                       exact = NA, correct = TRUE) %>% 
  rownames_to_column() %>% 
  select(rowname, pvalue)

#Categorical
chisq.test(table(meta_FULL$tb, meta_FULL$sex)) %>% tidy()
chisq.test(table(meta_FULL$tb, meta_FULL$race)) %>% tidy()

chisq.test(table(meta_FULL$tb, meta_FULL$smokhx2)) %>% tidy()
chisq.test(table(meta_FULL$tb, meta_FULL$bl_hiv)) %>% tidy()
chisq.test(table(meta_FULL$tb, meta_FULL$diabetes)) %>% tidy()

#### Statistics: case-control cohort ####
wilcox.x2 <- meta_CC %>% 
  filter(tb == "contact") %>% 
  select(bl_age, bmi) 
wilcox.y2 <- meta_CC %>% 
  filter(tb == "culture_confirm") %>% 
  select(bl_age, bmi) 

#Numeric
col_wilcoxon_twosample(wilcox.x2, wilcox.y2, alternative = "two.sided",
                       exact = NA, correct = TRUE) %>% 
  rownames_to_column() %>% 
  select(rowname, pvalue)

#Categorical
chisq.test(table(meta_CC$tb, meta_CC$sex)) %>% tidy()
chisq.test(table(meta_CC$tb, meta_CC$race)) %>% tidy()
chisq.test(table(meta_CC$tb, meta_CC$smokhx2)) %>% tidy()
# chisq.test(table(meta_CC$tb, meta_CC$bl_hiv)) %>% tidy()
# chisq.test(table(meta_CC$tb, meta_CC$diabetes)) %>% tidy()


#### A better case-control ####
meta_CC2 <- meta_CC %>% 
  # filter(smokhx2=="never") %>%
  filter(bl_age>=10) %>% 
  filter(!is.na(diabetes))

meta_CC %>% count(tb)
meta_CC2 %>% count(tb)
meta_CC2 %>% count(tb, smokhx2)

wilcox.x2 <- meta_CC2 %>% 
  filter(tb == "contact") %>% 
  select(bl_age, bmi) 
wilcox.y2 <- meta_CC2 %>% 
  filter(tb == "culture_confirm") %>% 
  select(bl_age, bmi) 

#Numeric
col_wilcoxon_twosample(wilcox.x2, wilcox.y2, alternative = "two.sided",
                       exact = NA, correct = TRUE) %>% 
  rownames_to_column() %>% 
  select(rowname, pvalue)

#Categorical
chisq.test(table(meta_CC2$tb, meta_CC2$sex)) %>% tidy()
chisq.test(table(meta_CC2$tb, meta_CC2$race)) %>% tidy()
