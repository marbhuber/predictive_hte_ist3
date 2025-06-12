##############################################################################
# Predictive HTE
# IST-3 trial
# 12 June 2025
# Markus Huber (markus.huber@insel.ch)
# The code is heavily influenced by Ewout Steyerberg's blog entry
# "Implementation of the PATH Statement" (https://www.fharrell.com/post/path/)
#  However, we follow the steps outline in Table 2 of the PATH statement by
# using a "treatment blinded" model for baseline risk
##############################################################################

##########################
# packages and environment

rm(list=ls())
library(tidyverse)
library(compareGroups)
library(labelled)
library(flextable)
library(emmeans)
library(rms)
library(DescTools)
library(cowplot)
library(haven)
library(ggridges)
options(scipen=999)

###########
# read data
# note that the data can be downloaded from the University of Edinburgh Datashare server
# https://datashare.ed.ac.uk/handle/10283/1931

data <- read_sas("C:/My Programs KAS/KAS Research/M. Huber/HTE SHORT/data/datashare_aug2015.sas7bdat")

df <- data %>% 
  mutate(
    gcs_eye_rand    = factor(gcs_eye_rand),
    gcs_motor_rand  = factor(gcs_motor_rand),
    gcs_verbal_rand = factor(gcs_verbal_rand)
  ) %>% 
  transmute(
    pid               = randhosp_id,
    treatment.factor  = factor(treatment,levels = c("Placebo","rt-PA")),
    treatment.numeric = as.integer(factor(treatment,levels = c("Placebo","rt-PA")))-1,
    age,
    sex               = factor(gender,levels = c(1,2),labels = c("Female","Male")),
    livealone_rand    = factor(livealone_rand,levels = c(1,2),labels = c("Yes","No")),
    infarct           = factor(infarct,levels = c(0,1,2),labels = c('No','Possibly Yes','Definitely Yes')),
    antiplat_rand     = factor(antiplat_rand,levels = c(1,2),labels = c("Yes","No")),
    atrialfib_rand    = factor(atrialfib_rand,levels = c(1,2),labels = c("Yes","No")),
    sbprand,
    dbprand,
    weight,
    gcseye = factor(
      case_when(
        gcs_eye_rand==1~"Non-spontaneously", # Never
        gcs_eye_rand==2~"Non-spontaneously", # To pain
        gcs_eye_rand==3~"Non-spontaneously", # To command
        gcs_eye_rand==4~"Spontaneously",
        TRUE~NA_character_)),
    gcsmotor = factor(
      case_when(
        gcs_motor_rand==1~"Not normal", # None
        gcs_motor_rand==2~"Not normal", # Extend to pain
        gcs_motor_rand==3~"Not normal", # Abnormal flex to pain
        gcs_motor_rand==4~"Not normal", # normal flex to pain
        gcs_motor_rand==5~"Not normal", # Localises movements to pain
        gcs_motor_rand==6~"Normal",
        TRUE~NA_character_)),
    gcsverbal = factor(
      gcs_verbal_rand,levels = c(1:5),
      labels = c('None','Noises only','Inappropriate words','Confused in time, place or person','Orientated in time, place and person')
    ),
    gcs_score_rand,
    nihss,
    stroketype = factor(case_when(
      stroketype==1~"TACI",
      stroketype==2~"PACI",
      stroketype==3~"LACI",
      stroketype==4~"POCI",
      TRUE~NA)), # Note that this excluded "Other" stroketypes (only 5 patients)
    # Primary outcome: Alive and independent at 6 months (OHS 0-2)
    outcome   = case_when( 
      aliveind6==1~1, # value yn: 1='Yes'
      aliveind6==2~0, # value yn: 2='No'
    )) %>%
  na.omit() %>% # complete case analysis
  # define subgroups
  mutate(
    subgroup.age            = case_when(age<=80~"Subgroup 1",age>80~"Subgroup 2",TRUE~NA),
    subgroup.sex            = case_when(sex=="Female"~"Subgroup 1",sex=="Male"~"Subgroup 2",TRUE~NA))

################
# Export Table 1
################

# summary measures
export2word(
createTable(
  compareGroups(~.-pid,data=df,method=NA,simplify = F),show.n=T,digits = 1
),file = "table1_raw_12june2025.docx"
)

# confidence intervals (for primary outcome)
export2word(
  createTable(
    compareGroups(~.-pid,data=df,method=NA,simplify = F),show.n=T,digits = 1,show.ci=T
  ),file = "table1_raw_12june2025_ci.docx"
)

################
################
# Risk modelling
################
################

# logistic regression model for Table 1
# it does not includes treatment variable as in Table 2 of the PATH statement

mymod <- glm(outcome ~ 
               age+
               sex+
               livealone_rand+
               infarct+
               antiplat_rand+
               atrialfib_rand+
               sbprand+
               dbprand+
               weight+
               gcseye+
               gcsmotor+
               gcsverbal+
               gcs_score_rand+
               nihss+
               stroketype,
             data=df,
             family=binomial())

save_as_docx(gtsummary::tbl_regression(mymod,exponentiate=T) %>% 
  gtsummary::as_flex_table(),
  path = "regression_model_raw_12june25.docx")

####################################################################
# From now - following Steyerberg's analysis - use this fitted model
# without treatment term - "blinded"

risk.mod <- lrm(outcome ~ 
                 age+
                 sex+
                 livealone_rand+
                 infarct+
                 antiplat_rand+
                 atrialfib_rand+
                 sbprand+
                 dbprand+
                 weight+
                 gcseye+
                 gcsmotor+
                 gcsverbal+
                 gcs_score_rand+
                 nihss+
                 stroketype,
               data=df,
               x=T,
               maxit=99)

# linear predictor of baseline risk and attach to data frame
linear.predictor.baseline <- risk.mod$linear.predictors
# linear predictor on log odds scale
df$linear.predictor.baseline  <- as.vector(linear.predictor.baseline)
# calculate baseline risks on probability scale
df$predictions                <- plogis(as.vector(linear.predictor.baseline))

#########################################################
#########################################################
# Evaluate interaction of treatment with linear predictor

mod.interaction <- lrm(outcome ~ treatment.factor * linear.predictor.baseline, data=df)
print(anova(mod.interaction)) # p = 0.0003

#############
#############
# Risk groups

groups              <- cut2(linear.predictor.baseline, g=4)
df$groups.quartiles <- groups

df <- df %>% 
  mutate(
    groups.quartiles.numeric = as.numeric(groups.quartiles)
  ) %>% 
  mutate(
    subgroup.risk = case_when(
      groups.quartiles.numeric==1~"Subgroup 1", # lowest risk
      groups.quartiles.numeric==2~"Subgroup 2",
      groups.quartiles.numeric==3~"Subgroup 3",
      groups.quartiles.numeric==4~"Subgroup 4"  # highest risk
    )
  )

####################################
####################################
# Treatment benefits in risk groups
# group0:   control group (Placebo)
# group1:   treatment group (rt-PA)
# ratediff: treatment minus control
# since we have a beneficial outcome
#   ratediff>0 favors treatment
#   ratediff<0 favors control

df.benefits <- data.frame()

groups    <- cut2(linear.predictor.baseline, g=4)
group0    <- groups[df$treatment.numeric==0]
group1    <- groups[df$treatment.numeric==1]
rate0     <- prop.table(table(group0, df$outcome[df$treatment.numeric==0]),1 )[,2]
rate1     <- prop.table(table(group1, df$outcome[df$treatment.numeric==1]),1 )[,2]
ratediff  <- rate1-rate0
events0   <- table(group0, df$outcome[df$treatment.numeric==0])[,2] # events in control group (placebo)
nevents0  <- table(group0, df$outcome[df$treatment.numeric==0])[,1] # no events in control group (placebo)
events1   <- table(group1, df$outcome[df$treatment.numeric==1])[,2] # events in treatment group (rt-PA)
nevents1  <- table(group1, df$outcome[df$treatment.numeric==1])[,1] # no events in treatment group (rt-PA)
n0        <- events0 + nevents0
n1        <- events1 + nevents1
CI        <- BinomDiffCI(events1, n1, events0, n0, method = "scorecc")
mean.risk <- df %>% group_by(groups.quartiles) %>% summarise(mean.risk = mean(predictions)) %>% pull(mean.risk) %>% as.numeric()
df.benefits <- rbind(df.benefits,data.frame(
  Subgroup  = c("Q1","Q2","Q3","Q4"),
  risk      = mean.risk,
  mymean    = ratediff,
  mylower   = as.numeric(CI[,2]),
  myupper   = as.numeric(CI[,3])
))

# age
groups    <- df$subgroup.age
group0    <- groups[df$treatment.numeric==0]
group1    <- groups[df$treatment.numeric==1]
rate0     <- prop.table(table(group0, df$outcome[df$treatment.numeric==0]),1 )[,2]
rate1     <- prop.table(table(group1, df$outcome[df$treatment.numeric==1]),1 )[,2]
ratediff  <- rate1-rate0
events0   <- table(group0, df$outcome[df$treatment.numeric==0])[,2]
nevents0  <- table(group0, df$outcome[df$treatment.numeric==0])[,1]
events1   <- table(group1, df$outcome[df$treatment.numeric==1])[,2]
nevents1  <- table(group1, df$outcome[df$treatment.numeric==1])[,1]
n0        <- events0 + nevents0
n1        <- events1 + nevents1
CI        <- BinomDiffCI(events1, n1, events0, n0, method = "scorecc")
mean.risk <- df %>% group_by(subgroup.age) %>% summarise(mean.risk = mean(predictions)) %>% pull(mean.risk) %>% as.numeric()
df.benefits <- rbind(df.benefits,data.frame(
  Subgroup   = c("\U2264 80 yrs","> 80 yrs"), # Subgroup 1: below 80, Subgroup 2: above 80
  risk    = mean.risk,
  mymean  = ratediff,
  mylower = as.numeric(CI[,2]),
  myupper = as.numeric(CI[,3])
))

# sex
groups    <- df$subgroup.sex
group0    <- groups[df$treatment.numeric==0]
group1    <- groups[df$treatment.numeric==1]
rate0     <- prop.table(table(group0, df$outcome[df$treatment.numeric==0]),1 )[,2]
rate1     <- prop.table(table(group1, df$outcome[df$treatment.numeric==1]),1 )[,2]
ratediff  <- rate1-rate0
events0   <- table(group0, df$outcome[df$treatment.numeric==0])[,2]
nevents0  <- table(group0, df$outcome[df$treatment.numeric==0])[,1]
events1   <- table(group1, df$outcome[df$treatment.numeric==1])[,2]
nevents1  <- table(group1, df$outcome[df$treatment.numeric==1])[,1]
n0        <- events0 + nevents0
n1        <- events1 + nevents1
CI        <- BinomDiffCI(events1, n1, events0, n0, method = "scorecc")
mean.risk <- df %>% group_by(subgroup.sex) %>% summarise(mean.risk = mean(predictions)) %>% pull(mean.risk) %>% as.numeric()
df.benefits <- rbind(df.benefits,data.frame(
  Subgroup   = c("Female","Male"), # Subgroup 1: Female, Subgroup 2: Male
  risk    = mean.risk,
  mymean  = ratediff,
  mylower = as.numeric(CI[,2]),
  myupper = as.numeric(CI[,3])
)) %>% 
  mutate(
    Subgroup = factor(Subgroup,levels = c("\U2264 80 yrs","> 80 yrs","Female","Male","Q1","Q2","Q3","Q4")))

df.benefits <- df.benefits %>% 
  mutate(
    type = case_when(
      Subgroup == "\U2264 80 yrs" ~ "one-variable-at-a-time subgroups",
      Subgroup == "> 80 yrs"~ "one-variable-at-a-time subgroups",
      Subgroup == "Male"~ "one-variable-at-a-time subgroups",
      Subgroup == "Female"~ "one-variable-at-a-time subgroups",
      Subgroup == "Q1"~"risk-based subgroups",
      Subgroup == "Q2"~"risk-based subgroups",
      Subgroup == "Q3"~"risk-based subgroups",
      Subgroup == "Q4"~"risk-based subgroups"
    )
  )

##########################
##########################
# density of baseline risk

df.density    = data.frame()
xp            = seq(0.002,0.9,by=0.001) # up 90% baseline risk
density       = density(df$predictions,n=1e3)
df.density    <- rbind(
  df.density,
  data.frame(
    risk    = density$x,
    mylower = -0.12,
    myupper = -0.12+(20e-3*density$y)
  ) %>% 
    filter(risk<0.9)
)

###################################################
# Average treatment effect with logistic regression
# Contrast: rt-PA minus Placebo
# Such that contrast is consistent with ratediff

mod.ate <- glm(outcome~treatment.factor,data = df,family = binomial())
df.ate  <- emmeans(
  mod.ate,~treatment.factor,transform = "response") %>% 
  pairs(infer = c(T,T),reverse=T) %>% 
  as.data.frame() %>% 
  mutate(
    mylabel = paste0(
      "Average treatment effect:\n",
      format(round(100*estimate,1),nsmall = 1),
      "%\n(95%-CI: ",
      format(round(100*asymp.LCL,1),nsmall = 1),
      "% to ",
      format(round(100*asymp.UCL,1),nsmall = 1),
      "%)"
    )
  )

######################
# median and mean risk
mean.risk        = df$predictions %>% mean()
median.risk      = df$predictions %>% median()
index.meanrisk   = which.min(abs(xp - mean.risk))
index.medianrisk = which.min(abs(xp - median.risk))

########
########
########
# Figure

figA <- df %>% 
  pivot_longer(c(subgroup.risk,subgroup.age,subgroup.sex)) %>% 
  mutate(
    ylab = factor(case_when(
      name=="subgroup.age"~"Age",
      name=="subgroup.sex"~"Sex",
      name=="subgroup.risk" ~ "Risk quartiles"),
    levels = c(
      "Age",
      "Sex",
      "Risk quartiles"))) %>% 
  select(ylab,value,predictions) %>% 
  mutate(
    Subgroup = factor(case_when(
      (ylab=="Risk quartiles" & value=="Subgroup 1")~"Q1",
      (ylab=="Risk quartiles" & value=="Subgroup 2")~"Q2",
      (ylab=="Risk quartiles" & value=="Subgroup 3")~"Q3",
      (ylab=="Risk quartiles" & value=="Subgroup 4")~"Q4",
      (ylab=="Age" & value=="Subgroup 1")~"\U2264 80 yrs",
      (ylab=="Age" & value=="Subgroup 2")~"> 80 yrs",
      (ylab=="Sex" & value=="Subgroup 1")~"Male",
      (ylab=="Sex" & value=="Subgroup 2")~"Female"),
    levels = c("\U2264 80 yrs","> 80 yrs","Male","Female","Q1","Q2","Q3","Q4"))) %>% 
  ggplot(aes(y = ylab)) +
  geom_density_ridges(
    aes(x = predictions, fill = Subgroup),
    alpha = .4, color = "white", height=0.3
  )+
  scale_x_continuous(labels = scales::percent,n.breaks = 10,limits = c(0,0.8))+
  coord_cartesian(clip = "off") +
  # theme_ridges(grid = FALSE)+
  theme_classic()+
  ggsci::scale_fill_aaas()+#nejm()+
  ylab("")+
  xlab("Baseline probability of being alive and independent at 6-months follow-up")+
  theme(
    legend.position = "bottom",legend.direction = "horizontal", #c(0.9,0.15),
    text = element_text(size=14),
    legend.key.size = unit(.5, 'cm'), #change legend key size
    legend.key.height = unit(.5, 'cm'), #change legend key height
    legend.key.width = unit(.5, 'cm'), #change legend key width
    legend.title = element_text(size=11), #change legend title font size
    legend.text = element_text(size=10))+
  guides(fill = guide_legend(nrow = 1))+
  labs(subtitle = "Risk distribution in selected subgroups")

figB <- ggplot()+
  geom_rect(data=df.ate, mapping=aes(xmin=-Inf, xmax=Inf, ymin=asymp.LCL, ymax=asymp.UCL), fill="black",color=NA, alpha=0.05)+
  geom_hline(yintercept = df.ate$estimate,color="black",linetype = "dashed")+
  annotate('text', x = 0.85, y = df.ate$estimate+0.007,label = 'ATE', 
           size = 3.7)+
  geom_ribbon(data=df.density,aes(x=risk,ymin=mylower,ymax=myupper),colour=NA,fill="grey",show.legend = F,alpha=0.7)+
  geom_point(data=df.benefits,aes(x=risk,y=mymean,color=Subgroup),size=4)+
  geom_errorbar(data = df.benefits,aes(x=risk,y=mymean,ymin=mylower,ymax=myupper,color=Subgroup),width=0.012,alpha=0.6)+
  geom_text(data=df.benefits,aes(x=risk,y=mymean,label=Subgroup,color=Subgroup),vjust = -1.2,hjust=-.1,show.legend = F)+
  geom_point(data = data.frame(risk=xp[index.meanrisk]),aes(x=risk,y=-0.12),shape = 2,size=2.5,fill="black")+
  annotate(geom="text", x=0.46, y=-0.12, label="average risk",color="black",size=3)+
  geom_point(data = data.frame(risk=xp[index.medianrisk]),aes(x=risk,y=-0.12),shape = 0,size=2.5,fill="black")+
  annotate(geom="text", x=0.15, y=-0.12, label="median (typical) risk",color="black",size=3)+
  annotate('text', x = 0, y = 0.07,  
           label = 'benefit', 
           size = 3.7, 
           angle='90')+
  annotate('text', x = 0, y = -0.07,  
           label = 'harm', 
           size = 3.7, 
           angle='90')+
  scale_y_continuous(labels = scales::percent,breaks = seq(-.12,.14,by=.02))+
  scale_x_continuous(labels = scales::percent,breaks = seq(0,0.9,by=0.1),limits = c(0,0.9))+
  theme_classic()+
  theme(
    text = element_text(size=14),
    legend.title = element_blank(),
    legend.position="none",#bottom",
    legend.direction="horizontal")+
  xlab("Baseline probability of being alive and independent at 6-months follow-up")+
  ylab("Treatment benefit\n(rt-PA versus Placebo)")+
  guides(color = guide_legend(nrow = 1))+
  ggsci::scale_color_aaas()+
  facet_wrap(~type)

cowplot::plot_grid(
  figA,
  figB,
  nrow=2,
  labels = "AUTO",
  rel_heights = c(1.75,2),
  align = "hv"
)

ggsave("figure1_12june2025.jpeg",dpi=600,width=8,height=8)
