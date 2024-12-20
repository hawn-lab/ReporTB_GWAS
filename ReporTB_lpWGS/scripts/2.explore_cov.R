library(tidyverse)
library(broom)
library(patchwork)

#### Data ####
load("~/project/data/metadata/ReportTB_metadata_clean.RData")

#### Summary full ####
meta_FULL %>% count(tb)

# Follow-up
meta_FULL %>% 
  group_by(tb) %>% 
  summarise(mean=mean(days_followup, na.rm=TRUE),
            sd=sd(days_followup, na.rm=TRUE))
t.test(days_followup~tb, data=meta_FULL) %>% tidy()

# Age
meta_FULL %>% 
  group_by(tb) %>% 
  summarise(mean=mean(bl_age, na.rm=TRUE),
            sd=sd(bl_age, na.rm=TRUE))
t.test(bl_age~tb, data=meta_FULL) %>% tidy()

# Sex
meta_FULL %>% count(tb, sex) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_FULL$tb, meta_FULL$sex)) %>% tidy()

# Race
meta_FULL %>% count(tb, race) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_FULL$tb, meta_FULL$race)) %>% tidy()

# BMI
meta_FULL %>% 
  group_by(tb) %>% 
  summarise(mean=mean(bmi, na.rm=TRUE),
            sd=sd(bmi, na.rm=TRUE))
t.test(bmi~tb, data=meta_FULL) %>% tidy()

# Smoking
# meta_FULL %>% count(tb, smokhx) %>%
#   group_by(tb) %>% mutate(perc = n/sum(n))
# chisq.test(table(meta_FULL$tb, meta_FULL$smokhx)) %>% tidy()

meta_FULL  %>% count(tb, smokhx2) %>%
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_FULL$tb, meta_FULL$smokhx2)) %>% tidy()

# HIV
meta_FULL %>% count(tb, bl_hiv) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_FULL$tb, meta_FULL$bl_hiv)) %>% tidy()

# Diabetes
meta_FULL %>% count(tb, diabetes) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_FULL$tb, meta_FULL$diabetes)) %>% tidy()

# QTF
meta_FULL  %>% 
  count(tb, qtf) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_FULL$tb, meta_FULL$qtf)) %>% tidy()

#### Summary CC ####
meta_CC %>% count(tb)

# Follow-up
meta_CC %>% 
  group_by(tb) %>% 
  summarise(mean=mean(days_followup, na.rm=TRUE),
            sd=sd(days_followup, na.rm=TRUE))
t.test(days_followup~tb, data=meta_CC) %>% tidy()

# Age
meta_CC %>% 
  group_by(tb) %>% 
  summarise(mean=mean(bl_age, na.rm=TRUE),
            sd=sd(bl_age, na.rm=TRUE))
t.test(bl_age~tb, data=meta_CC) %>% tidy()

# Sex
meta_CC %>% count(tb, sex) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_CC$tb, meta_CC$sex)) %>% tidy()

# Race
meta_CC %>% count(tb, race) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_CC$tb, meta_CC$race)) %>% tidy()

# BMI
meta_CC %>% 
  group_by(tb) %>% 
  summarise(mean=mean(bmi, na.rm=TRUE),
            sd=sd(bmi, na.rm=TRUE))
t.test(bmi~tb, data=meta_CC) %>% tidy()

# Smoking
# meta_CC %>% count(tb, smokhx) %>% 
#   group_by(tb) %>% mutate(perc = n/sum(n))
# chisq.test(table(meta_CC$tb, meta_CC$smokhx)) %>% tidy()

meta_CC %>% count(tb, smokhx2) %>%
  group_by(tb) %>% mutate(perc = n/sum(n))
chisq.test(table(meta_CC$tb, meta_CC$smokhx2)) %>% tidy()

# HIV
meta_CC %>% count(tb, bl_hiv) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
# chisq.test(table(meta_CC$tb, meta_CC$bl_hiv)) %>% tidy()

# Diabetes
meta_CC %>% count(tb, diabetes) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
# chisq.test(table(meta_CC$tb, meta_CC$diabetes)) %>% tidy()

# QTF
meta_CC  %>% 
  count(tb, qtf) %>% 
  group_by(tb) %>% mutate(perc = n/sum(n))
# chisq.test(table(meta_CC$tb, meta_CC$qtf)) %>% tidy()


#### Plots ####
#Numeric
p1 <- meta_FULL %>% 
  select(sample_id, tb,bl_age, bmi, days_followup) %>% 
  mutate(tb=recode(tb, "culture_confirm"="tb")) %>% 
  pivot_longer(-c(sample_id,tb)) %>% 
  
  ggplot(aes(x=tb,y=value)) +
  geom_jitter(width=0.2, height=0, color="grey") +
  #Stdev bar
  stat_summary(fun.data=mean_sdl, 
               fun.args = list(mult=1), 
               geom="errorbar", width=0.25) +
  #mean bar
  stat_summary(fun=mean, geom="errorbar",
               aes(ymax=after_stat(y), ymin=after_stat(y)),
               width=0.5) + 
  theme_classic() +
  labs(x="") +
  facet_wrap(~name, nrow=1, scales="free")

#Categorical
cat_ls <- list()

for(var in c("sex", "race", "bl_hiv", "smokhx", "diabetes")){
  temp <- meta_FULL %>% 
    #make other race group
    mutate(race = case_when(race %in% c("asian","indian",NA,"6") ~ "other",
                            TRUE ~ race)) %>% 
    select(sample_id, tb, all_of(var)) %>% 
    pivot_longer(-c(sample_id,tb)) %>% 
    count(tb, name, value) %>% 
    ungroup() %>% 
    mutate(tb=recode(tb, "culture_confirm"="tb")) %>% 
    ggplot(aes(x=tb,y=n, fill=value)) +
    geom_bar(position = "fill", stat = "identity") +
    theme_classic() +
    labs(x="", y="Proportion of participants", fill="") +
    facet_wrap(~name, nrow=1, scales="free")
  
  #recolor sex
  if(var=="sex"){
    cat_ls[[var]] <- temp +
      scale_fill_manual(values=c("lightblue","orange"))
  } else {
    cat_ls[[var]] <- temp 
  }
}

lo <- "
AABC
DEFF"
plot_all <- p1 + cat_ls[c(2,4)] + cat_ls[c(1,3,5)] +
  plot_layout(design = lo)
# plot_all

ggsave("result/meta_plots.png", plot = plot_all, width = 10, height = 6)
