#################################################################################################################################
#### Data analysis for "Behavioural variation and plasticity along an invasive ant introduction pathway" (Felden et al. 2018) ###
#################################################################################################################################
#### Contact: antoine@felden.fr

#### Working directory, data files, libraries ###################################################################################

rm(list = ls())
DataDir = "/Users/antoinefelden/Documents/Research/Manuscripts/4-JAnimalEcol/your_directory/01_data"
FigDir = "/Users/antoinefelden/Documents/Research/Manuscripts/4-JAnimalEcol/your_directory/03_figures"

library(grid)
library(lme4)
library(multcomp)
library(survival)
library(ggplot2)
library(Hmisc)
library(Rmisc)
library(gridExtra)
library(ordinal)
library(pgirmess)

#### Colours ####################################################################################################################

black = rgb(0.5,0.5,0.5,1)
blue = rgb(0,0,0.5,1)
green = rgb(0,0.5,0,1)
red = rgb(0.5,0,0,1)

#### Data files #################################################################################################################

# Foraging activity
AR_file = read.csv(paste(DataDir,"/tree_ar.csv", sep = ""), h = T, sep = ",")
OZ_file = read.csv(paste(DataDir,"/tree_oz.csv", sep = ""), h = T, sep = ",")
CA_file = read.csv(paste(DataDir,"/tree_ca.csv", sep = ""), h = T, sep = ",")
NZ_file = read.csv(paste(DataDir,"/tree_nz.csv", sep = ""), h = T, sep = ",")
unfiltered_raw_data = rbind(AR_file, OZ_file,CA_file,NZ_file)

# Aggression
agg_det = read.csv(paste(DataDir,"/aggression_all.csv",sep=""), h = T, sep = ",")
agg_det <- na.omit(agg_det)

# HPLC
data_HPLC <- read.csv(paste(DataDir,"/Exp_1-pathway.csv",sep=""),header=T,na.strings=NA,
                      colClasses=c(rep("factor",5),rep("numeric",5),"factor"))

#### Variables ##################################################################################################################

unfiltered_raw_data$early_explo <- with(unfiltered_raw_data, (n_0 + n_10 + n_20))
unfiltered_raw_data$subcolony <- with(unfiltered_raw_data, paste(region, subcolony, sep = "")) # experimental colony ID
unfiltered_raw_data$id <- as.factor(seq(1:length(unfiltered_raw_data$day))) # id = a factor with one level per observation
unfiltered_raw_data$time_top2 <- ifelse(unfiltered_raw_data$time_top == 61,60,unfiltered_raw_data$time_top)
unfiltered_raw_data$status <- ifelse(unfiltered_raw_data$time_top == 61,0,1)

#### Functions ##################################################################################################################

# Function multiplot to plot multiple plots
# https://gist.github.com/sckott/8444444

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, 
                      labs=list(), labpos=list(c(0.5,0.03), c(0.03,0.5))) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
    
    if(!length(labs) == 0){
      grid.text(labs[1], x=labpos[[1]][1], y=labpos[[1]][2], gp=gpar(fontsize=16))
      grid.text(labs[2], x=labpos[[2]][1], y=labpos[[2]][2], rot=90, gp=gpar(fontsize=16))
    }
  }
}

# Get legend for multiple plotting
# https://stackoverflow.com/questions/12539348/ggplot-separate-legend-and-plot

library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# multple plot with ggplot (alternative)
# http://rstudio-pubs-static.s3.amazonaws.com/2852_379274d7c5734f979e106dcf019ec46c.html

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# Calculate effect size in %
calc_increase <- function(ref,x) {
  inc = x-ref
  perc_inc = inc/ref
  return(perc_inc)
}

#### Main data subsets ##########################################################################################################

# Regions
AR = droplevels(subset(unfiltered_raw_data, region == "AR"))
OZ = droplevels(subset(unfiltered_raw_data, region == "AU"))
CA = droplevels(subset(unfiltered_raw_data, region == "CA"))
NZ = droplevels(subset(unfiltered_raw_data, region == "NZ"))

# Phases
A_phase = subset(unfiltered_raw_data, day == "1" | day == "2" | day == "3" | day == "4"| day == "5")
T_phase = subset(unfiltered_raw_data, day == "8" | day == "9"| day == "10"| day == "11"| day == "12"| day == "13"| day == "14"| day == "15")
TR_phase = subset(unfiltered_raw_data, phase == "treward")

exp_raw_data = rbind(T_phase,TR_phase)

#################################################################################################################################
#### Statistical analysis #######################################################################################################
#################################################################################################################################

#### GLMMs treatment and region on foraging activity

## Exploration phase
glmer_T_int = glmer(early_explo ~ treatment * region + (1 | site/subcolony) + (1 | day), data = T_phase, family = poisson)
summary(glmer_T_int)

glmer_T = glmer(early_explo ~ treatment + region + (1 | site/subcolony) + (1 | day), data = T_phase, family = poisson)
summary(glmer_T)

lrt <- function (obj1, obj2) { 
  L0 <- logLik(obj1) 
  L1 <- logLik(obj2) 
  L01 <- as.vector(2 * abs(L0 - L1))
  df <- abs(attr(L1, "df") - attr(L0, "df"))
  list(L01 = L01, df = df, 
       "p-value" = pchisq(L01, df, lower.tail = FALSE)) 
} 

lrt(glmer_T_int, glmer_T ) 
lrt(glmer_T, glmer_T_int) 

glmer_T = glmer(early_explo ~ treatment * region + (1 | site/subcolony) + (1 | day), data = T_phase, family = poisson)
summary(glmer_T)

# effect size for AU
calc_increase(log(2.682846),log(2.682846-0.596759))
# Effect size for low_sugar in CA
calc_increase(log(2.682846+0.285919),log(2.682846+0.285919-0.588781))
# effect size for OA
calc_increase(log(2.682846),log(2.682846+0.799798))
# effect size for LSOA
calc_increase(log(2.682846),log(2.682846+1.282161))

## Exploitation phase
glmer_TR = glmer(early_explo ~ treatment * region + (1 | site/subcolony) + (1 | day), data = TR_phase, family = poisson)
summary(glmer_TR)

# effect size for AU
calc_increase(log(3.69991),log(3.69991-0.98529))
# effect size for low_sugar in AR (no significantly different elsewhere)
calc_increase(log(3.69991),log(3.69991+0.6791))
#effect size for OA in AR (no significant elsewhere but OZ)
calc_increase(log(3.69991),log(3.69991+0.8727))
# effect size for LSOA in AR (no significant elsewhere but OZ)
calc_increase(log(3.69991),log(3.69991+1.11481))


## Fligner-Killeen: Variance among regions

fligner.test(early_explo ~ region, data = subset(exp_raw_data, treatment == "control" & phase == "treward"))
ARCA=rbind(AR,CA)
AROZ=rbind(AR,OZ)
ARNZ=rbind(AR,NZ)
CAOZ=rbind(CA,OZ)
CANZ=rbind(CA,NZ)
OZNZ=rbind(OZ,NZ)
fligner.test(early_explo ~ region, data = subset(ARCA, treatment == "control" & phase == "treatment")) #NS
fligner.test(early_explo ~ region, data = subset(AROZ, treatment == "control" & phase == "treatment")) #p=0.091
fligner.test(early_explo ~ region, data = subset(ARNZ, treatment == "control" & phase == "treatment")) #NS
fligner.test(early_explo ~ region, data = subset(CAOZ, treatment == "control" & phase == "treatment")) #p<0.001
fligner.test(early_explo ~ region, data = subset(CANZ, treatment == "control" & phase == "treatment")) #p=0.047
fligner.test(early_explo ~ region, data = subset(OZNZ, treatment == "control" & phase == "treatment")) #p=0.087
fligner.test(early_explo ~ region, data = subset(ARCA, treatment == "control" & phase == "treward")) #NS
fligner.test(early_explo ~ region, data = subset(AROZ, treatment == "control" & phase == "treward")) #p=0.100
fligner.test(early_explo ~ region, data = subset(ARNZ, treatment == "control" & phase == "treward")) #NS
fligner.test(early_explo ~ region, data = subset(CAOZ, treatment == "control" & phase == "treward")) #p<0.01
fligner.test(early_explo ~ region, data = subset(CANZ, treatment == "control" & phase == "treward")) #NS
fligner.test(early_explo ~ region, data = subset(OZNZ, treatment == "control" & phase == "treward")) #p=0.079

## Time to the top of the tree (survival analysis)

testcox_all = coxph(Surv(time_top2, status) ~ treatment*region, data=T_phase)
print(summary(testcox_all))

# AR
survfit_obj_AR = survfit(Surv(time_top2, status) ~ treatment, data = subset(T_phase, region == "AR")) #LSOA: all colonies reached the top at 40' or before -> manual adjsutment (check survfit_obj$time)
survfit_coordinates_AR = data.frame(time=c(survfit_obj_AR$time,50,60),prop=c(1-survfit_obj_AR$surv,1,1),lower=c(1-survfit_obj_AR$upper,NA,NA),upper = c(1-survfit_obj_AR$lower,NA,NA))

# CA
survfit_obj_CA = survfit(Surv(time_top2, status) ~ treatment, data = subset(T_phase, region == "CA")) #not all colonies in any treatment reached the top at 60', no adjustment needed
survfit_coordinates_CA = data.frame(time=c(survfit_obj_CA$time),prop=c(1-survfit_obj_CA$surv),lower=c(1-survfit_obj_CA$upper),upper = c(1-survfit_obj_CA$lower))

# OZ
survfit_obj_OZ = survfit(Surv(time_top2, status) ~ treatment, data = subset(T_phase, region == "AU")) #LSOA: all colonies tothe top at 50' -> manual adjsutment
survfit_coordinates_OZ = data.frame(time=c(survfit_obj_OZ$time,60),prop=c(1-survfit_obj_OZ$surv,1),lower=c(1-survfit_obj_OZ$upper,NA),upper = c(1-survfit_obj_OZ$lower,NA))

# NZ
survfit_obj_NZ = survfit(Surv(time_top2, status) ~ treatment, data = subset(T_phase, region == "NZ")) #LSOA: all colonies to the top at 30' & OA: -> manual adjsutment
survfit_pre_coordinates_NZ = data.frame(time=c(survfit_obj_NZ$time,40,50,60),prop=c(1-survfit_obj_NZ$surv,1,1,1),lower=c(1-survfit_obj_NZ$upper,NA,NA,NA),upper = c(1-survfit_obj_NZ$lower,NA,NA,NA))
missing_line=survfit_pre_coordinates_NZ[15,]; missing_line$time <- 40
survfit_coordinates_NZ = rbind(survfit_pre_coordinates_NZ[1:15,],missing_line, survfit_pre_coordinates_NZ[16:23,])

# Bind all data and set up labels
survfit_coordinates_all = rbind(survfit_coordinates_AR,survfit_coordinates_CA,survfit_coordinates_OZ,survfit_coordinates_NZ)

treatment = c(rep("control",6),rep("low_sugar",6),rep("OA",6),rep("OA_low_sugar",6))
survfit_coordinates = data.frame(region = c(rep("AR",24),rep("CA",24),rep("AU",24),rep("NZ",24)),treatment = rep(treatment,4),survfit_coordinates_all)

zeros = data.frame(treatment = c("control","low_sugar","OA","OA_low_sugar"), time = c(0,0,0,0), prop = c(0,0,0,0),lower = c(0,0,0,0),upper = c(0,0,0,0))
zeros=rbind(zeros,zeros,zeros,zeros)
zeros_reg = cbind(region=c(rep("AR",4),rep("CA",4),rep("AU",4),rep("NZ",4)),zeros)

survfit_coordinates = rbind(zeros_reg,survfit_coordinates)
survfit_coordinates <- survfit_coordinates[order(survfit_coordinates$time),]
survfit_coordinates[is.na(survfit_coordinates)] <- 1

### Aggression (cumulative link model)

df_lines <- seq(1:nrow(agg_det))
agg_det$line <- df_lines
agg_det$subcolo <- interaction(agg_det$region,agg_det$subcolony)
agg_det$reac <- agg_det$reactive_attack+agg_det$chemicals
agg_det$n_agg = agg_det$proactive_attack+agg_det$reactive_attack+agg_det$chemicals
agg_det$n_av = agg_det$avoidance+agg_det$neutral
aggression_full = NULL
for (li in df_lines) {
  data_assay = subset(agg_det, line == li)
  df_assay=NULL
  if (data_assay$n_av != 0) {
    df_assay <- rbind(df_assay,data.frame("region" = rep(data_assay$region,data_assay$n_av), "subcolony" = rep(data_assay$subcolo,data_assay$n_av), "treatment" = rep(data_assay$treatment,data_assay$n_av), "phase" = rep(data_assay$phase,data_assay$n_av), "assay" = rep(data_assay$assay,data_assay$n_av), "agg_full" = rep(0,data_assay$n_av)))}
  if (data_assay$neutral != 0) {
    df_assay <- rbind(df_assay,data.frame("region" = rep(data_assay$region,data_assay$neutral), "subcolony" = rep(data_assay$subcolo,data_assay$neutral), "treatment" = rep(data_assay$treatment,data_assay$neutral), "phase" = rep(data_assay$phase,data_assay$neutral), "assay" = rep(data_assay$assay,data_assay$neutral), "agg_full" = rep(1,data_assay$neutral)))}
  if (data_assay$reac != 0) {
    df_assay <- rbind(df_assay,data.frame("region" = rep(data_assay$region,data_assay$reac), "subcolony" = rep(data_assay$subcolo,data_assay$reac), "treatment" = rep(data_assay$treatment,data_assay$reac), "phase" = rep(data_assay$phase,data_assay$reac), "assay" = rep(data_assay$assay,data_assay$reac), "agg_full" = rep(2,data_assay$reac)))}
  if (data_assay$proactive_attack != 0) {
    df_assay <- rbind(df_assay,data.frame("region" = rep(data_assay$region,data_assay$proactive_attack), "subcolony" = rep(data_assay$subcolo,data_assay$proactive_attack), "treatment" = rep(data_assay$treatment,data_assay$proactive_attack), "phase" = rep(data_assay$phase,data_assay$proactive_attack), "assay" = rep(data_assay$assay,data_assay$proactive_attack), "agg_full" = rep(3,data_assay$proactive_attack)))}
  if (data_assay$fight != 0) {
    df_assay <- rbind(df_assay,data.frame("region" = rep(data_assay$region,data_assay$fight), "subcolony" = rep(data_assay$subcolo,data_assay$fight), "treatment" = rep(data_assay$treatment,data_assay$fight), "phase" = rep(data_assay$phase,data_assay$fight), "assay" = rep(data_assay$assay,data_assay$fight), "agg_full" = rep(4,data_assay$fight)))}
  aggression_full= rbind(aggression_full,df_assay)}

clm_agg_rAR_int <- clm(as.factor(aggression_full$agg_full) ~ treatment * region,data= aggression_full, Hess=T,link = "logit", threshold = "flexible")
summary(clm_agg_rAR_int)

clm_agg_rAR <- clm(as.factor(aggression_full$agg_full) ~ treatment + region,data= aggression_full, Hess=T,link = "logit", threshold = "flexible")
summary(clm_agg_rAR)

### Biogenic amine titres (HPLC)

internal_stand = NULL
for (bat in unique(data_HPLC$Batch)) {
  sub_batch <- subset(data_HPLC, Batch == bat,select=c(Batch,DHBA,Buffer))
  line_int <- c(batch=bat,colony=sub_batch$Colony,mean_int=mean(na.omit(sub_batch$DHBA)),buffer=sub_batch$Buffer[1])
  internal_stand = rbind(internal_stand,line_int)
}
internal_stand <- data.frame(internal_stand)
internal_stand$mean_int <- as.integer(as.character(internal_stand$mean_int))

# Compute standard areas (samples and calibration)
concentrations <- data.frame(DA_SER = c(0,2.5,4,7,10), OA_TYR = c(0,1.35,2,3.5,5))

std_areas = NULL
for (buff in unique(data_HPLC$Buffer)) {
  standards <- subset(data_HPLC, Buffer == buff)
  
  for (bat in unique(standards$Batch)) {
    sub_batch <- cbind(subset(standards, Batch == bat))
    batch_index <- which(internal_stand$batch == bat)
    sub_batch$std_OA <- (internal_stand[batch_index,2]/sub_batch$DHBA)*sub_batch$OA
    sub_batch$std_DHBA <- (internal_stand[batch_index,2]/sub_batch$DHBA)*sub_batch$DHBA
    sub_batch$std_DA <- (internal_stand[batch_index,2]/sub_batch$DHBA)*sub_batch$DA
    sub_batch$std_TYR <- (internal_stand[batch_index,2]/sub_batch$DHBA)*sub_batch$TYR
    sub_batch$std_SER <- (internal_stand[batch_index,2]/sub_batch$DHBA)*sub_batch$SER
    std_areas <- rbind(std_areas,sub_batch)
  }
}

# Plot calibration curves and compute coefficients (in picograms/microliter - or picograms for amount values)
concentrations <- data.frame(conc_DA_SER = c(0,2.5,4,7,10), conc_OA_TYR = c(0,1.35,2,3.5,5), DA_SER = c(0,2.5*40,4*40,7*40,10*40), OA_TYR = c(0,1.35*40,2*40,3.5*40,5*40))

std_calibration <- subset(std_areas, Colony == "STD")
stds_amounts = NULL
for (bat in unique(std_calibration$Batch)) {
  sub_batch <- cbind(subset(std_calibration, Batch == bat),concentrations)
  stds_amounts <- rbind(stds_amounts,sub_batch)
}

calibration=NULL
for (buff in unique(stds_amounts$Buffer)) {
  sub_buff_stds <- subset(stds_amounts,Buffer == buff)
  par(mfrow=c(2,2))
  lm_OA <- lm(OA_TYR~std_OA,data=sub_buff_stds)
  plot(OA_TYR~std_OA,data=sub_buff_stds,main="OA")
  line_OA <- data.frame(a_OA = lm_OA$coefficients[1], b_OA = lm_OA$coefficients[2])
  abline(a=lm_OA$coefficients[1], b=lm_OA$coefficients[2])
  text(sub_buff_stds$std_OA,sub_buff_stds$OA_TYR,labels=sub_buff_stds$Batch,cex=2)
  lm_DA <- lm(DA_SER~std_DA,data=sub_buff_stds)
  plot(DA_SER~std_DA,data=sub_buff_stds,main="DA")
  line_DA <- data.frame(a_DA = lm_DA$coefficients[1], b_DA = lm_DA$coefficients[2])
  abline(a=lm_DA$coefficients[1], b=lm_DA$coefficients[2])
  text(sub_buff_stds$std_DA,sub_buff_stds$DA_SER,labels=sub_buff_stds$Batch,cex=2)
  lm_TYR <- lm(OA_TYR~std_TYR,data=sub_buff_stds)
  plot(OA_TYR~std_TYR,data=sub_buff_stds,main="TYR")
  line_TYR <- data.frame(a_TYR = lm_TYR$coefficients[1], b_TYR = lm_TYR$coefficients[2])
  abline(a=lm_TYR$coefficients[1], b=lm_TYR$coefficients[2])
  text(sub_buff_stds$std_TYR,sub_buff_stds$OA_TYR,labels=sub_buff_stds$Batch,cex=2)
  lm_SER <- lm(DA_SER~std_SER,data=sub_buff_stds)
  plot(DA_SER~std_SER,data=sub_buff_stds,main="SER")
  line_SER <- data.frame(a_SER = lm_SER$coefficients[1], b_SER = lm_SER$coefficients[2], buffer = buff)
  abline(a=lm_SER$coefficients[1], b=lm_SER$coefficients[2])
  text(sub_buff_stds$std_SER,sub_buff_stds$DA_SER,labels=sub_buff_stds$Batch,cex=2)
  calibration_buff <- data.frame(line_OA,line_DA,line_TYR,line_SER)
  calibration <- rbind(calibration, calibration_buff)
}

# Convert sample areas into compound amount (could have used predict.lm())
exp_samples <- subset(std_areas, Colony != "STD")
for (buff in unique(exp_samples$Buffer)) {
  buff_index <- which(calibration$buffer == buff)
  exp_samples$area_OA <- (calibration$a_OA[buff_index]+calibration$b_OA[buff_index]*exp_samples$std_OA)*((60/40)/5)
  exp_samples$area_DA <- (calibration$a_DA[buff_index]+calibration$b_DA[buff_index]*exp_samples$std_DA)*((60/40)/5)
  exp_samples$area_TYR <- (calibration$a_TYR[buff_index]+calibration$b_TYR[buff_index]*exp_samples$std_TYR)*((60/40)/5)
  exp_samples$area_SER <- (calibration$a_SER[buff_index]+calibration$b_SER[buff_index]*exp_samples$std_SER)*((60/40)/5)
}

exp_samples$amount_OA <- exp_samples$area_OA*((60/40)/5)
exp_samples$amount_DA <- exp_samples$area_DA*((60/40)/5)
exp_samples$amount_TYR <- exp_samples$area_TYR*((60/40)/5)
exp_samples$amount_SER <- exp_samples$area_SER*((60/40)/5)
exp_samples$OA_SER <- exp_samples$amount_OA/exp_samples$amount_SER


kruskal.test(amount_OA~Treatment,data=droplevels(exp_samples))
kruskalmc(amount_OA~Treatment,data=droplevels(exp_samples))

kruskal.test(amount_TYR~Treatment,data=droplevels(exp_samples))
kruskalmc(amount_TYR~Treatment,data=droplevels(exp_samples))

kruskal.test(amount_SER~Treatment,data=droplevels(exp_samples))
kruskalmc(amount_SER~Treatment,data=droplevels(exp_samples))

### Resource discovery (GLMM)

glmer_TR_sug = glmer(sug_10 ~ treatment + region  + (1 | site/subcolony) + (1 | day), data = TR_phase, family = poisson)
summary(glmer_TR_sug)

# Effect size for NZ
calc_increase(log(1.3865),log(1.3865+0.4768))
# Effect size for low_sugar
calc_increase(log(1.3865),log(1.3865+0.8526))
#effect size for OA
calc_increase(log(1.3865),log(1.3865+1.2615))
# effect size for LSOA
calc_increase(log(1.3865),log(1.3865+1.8249))

#################################################################################################################################
#### Figures ####################################################################################################################
#################################################################################################################################

#### THE FOLLOWING SECTION NEEDS TO BE RAN IN ORDER TO CORRECTLY PLOT FIGURES

### Data manipulation for plotting (relevel of factors)
exp_raw_data$region <- droplevels(factor(exp_raw_data$region,levels = c("AR","CA","AU","NZ"),ordered=T))
exp_raw_data$site <- factor(exp_raw_data$site,
                                   levels = c("OTA","ROS","UBA","VLO","BLB","CAL","OAK","RCH","CLI","MON","PHI","STA","HAS","PAR","PET","TAR"), ordered=T)

T_phase = subset(exp_raw_data, day == "8" | day == "9"| day == "10"| day == "11"| day == "12"| day == "13"| day == "14"| day == "15")
TR_phase = subset(exp_raw_data, phase == "treward")

### FIGURE 2: Variation in foraging in control groups along the introduction pathway

## Control groups, acclimation phase (not shown)
A_phase$region <- droplevels(factor(A_phase$region,levels = c("AR","CA","AU","NZ"),ordered=T))
ctrl_bplot_accli = ggplot(aes(y=early_explo,x=region),data=subset(A_phase, treatment == "control")) + theme_bw() +
  geom_boxplot(fill=c("cornflowerblue","tomato4","coral3","red1"),position=position_dodge(0.8),notch=T) +
  labs(title="0. Baseline Exploration of a novel environment",x="",y="Foraging activity \n (in number of workers foraging)\n") + ylim(0,180) +
  theme(axis.text.x = element_text(size = 14, colour = "black"),axis.text.y = element_text(size = 14, colour = "black"),axis.title.y = element_text(size = 16)) +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  scale_x_discrete(labels=c("Argentina","California","Australia","New Zealand")) +
  annotate("text", x = c(1,2,3,4), y=100, label = c("a", "a", "b", "a"),fontface = 2,size=6)
ctrl_bplot_accli

## Control groups, exploration phase
ctrl_bplot_explo = ggplot(aes(y=early_explo,x=region),data=subset(T_phase, treatment == "control")) + theme_bw() +
  geom_boxplot(fill=c("cornflowerblue","tomato4","coral3","red1"),position=position_dodge(0.8),notch=F) +
  labs(title="a. Exploration of a novel environment",x="",y="Foraging activity \n (in number of workers foraging)\n") + ylim(0,70) +
  theme(axis.text.x = element_text(size = 14, colour = "black"),axis.text.y = element_text(size = 14, colour = "black"),axis.title.y = element_text(size = 16)) +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  scale_x_discrete(labels=c("Argentina","California","Australia","New Zealand")) +
  annotate("text", x = c(1,2,3,4), y=70, label = c("a", "a", "b", "a"),fontface = 2,size=6)
ctrl_bplot_explo

## Control groups, exploitation phase
ctrl_bplot_exploit = ggplot(aes(y=early_explo,x=region),data=subset(TR_phase, treatment == "control")) + theme_bw() +
  geom_boxplot(fill=c("cornflowerblue","tomato4","coral3","red1"),position=position_dodge(0.8),notch=F) +
  labs(title="b. Exploitation of a new food source",x="",y="") + ylim(0,170) +
  theme(axis.text.x = element_text(size = 14, colour = "black"),axis.text.y = element_text(size = 14, colour = "black"),axis.title.y = element_text(size = 16)) +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  scale_x_discrete(labels=c("Argentina","California","Australia","New Zealand")) +
  annotate("text", x = c(1,2,3,4), y=170, label = c("a", "a", "b", "a"),fontface = 2,size=6)
ctrl_bplot_exploit

## Control groups broken down by sites, exploitation phase
sub_site_TRp <- data.frame(subset(exp_raw_data, phase == "treward" & treatment == "control"))
sub_site_TRp <- sub_site_TRp[order(sub_site_TRp$site),]
sub_site_TRp$site_replicates <- as.factor(rep(c(1,2),nrow(sub_site_TRp)/2))
site_bplot_TRp = ggplot(aes(y=early_explo,x=site,fill=region),data= sub_site_TRp) + theme_bw() +
  scale_fill_manual(values=c("cornflowerblue","tomato4","coral3","red1"),name="Region") + #,guide=F
  scale_color_manual(values=c("pink","cyan3"),name="Site \nreplication") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  labs(title="c. Exploitation of a new food source: variation among collection sites",x="",y="Foraging activity \n (in number of workers foraging)\n") +
  theme(plot.title=element_text(size=18, vjust=2),legend.position="right", legend.text=element_text(size=14),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.y=element_text(size = 16, colour = "black",vjust=1),
        axis.title.x=element_text(size = 16, colour = "black")) +
  geom_boxplot(outlier.size = 4, notch=F) +
  geom_point(aes(color = site_replicates),position = position_jitter(width = 0.2))+
  annotate("text",
           x = c(2.5, 6.5, 10.5, 14.5),
           y = rep(140,4),
           label = c("Argentina","California","Australia","New\nZealand"),family="",fontface = 4,size=4,
           family = "", fontface = 3, size=4)
site_bplot_TRp

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
print(ctrl_bplot_explo, vp = vplayout(1, 1))  # the big plot covers rows 1:2 and cols 1:3
print(ctrl_bplot_exploit, vp = vplayout(1, 2))
print(site_bplot_TRp, vp = vplayout(2, 1:2))

dev.copy2pdf(file=paste(FigDir,"/fig_panel.pdf",sep=""), 
             width=12, 
             height=8)

### FIGURE 3: Effect of treatments along the introduction pathway

## Exploration phase
T_bplot_reg = ggplot(aes(y=early_explo,x=region,fill=treatment),data=subset(exp_raw_data, phase == "treatment")) + theme_bw() +
  geom_boxplot(notch=F,position=position_dodge(0.88)) + ylim(0,190) + 
  scale_fill_manual(name="Treatment", labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="a. Exploration of a novel environment",x="",y="Foraging activity \n (in number of workers foraging)\n") +
  scale_x_discrete(labels=c("Argentina","California","Australia", "New Zealand")) +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="bottom", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,0,10),"pt")) +
  annotate("text", x = c(0.906,1.093,1.28,0.9995,1.9995,2.9995,3.9995),
         y=c(rep(140,3),rep(190,4)),
         label = c("NS","***","***","a","b","a","a"),fontface = 2,size=6) +
  geom_segment(aes(x=0.709,xend=1.28,y=175,yend=175)) +
  geom_segment(aes(x=1.709,xend=2.28,y=175,yend=175)) +
  geom_segment(aes(x=2.709,xend=3.28,y=175,yend=175)) +
  geom_segment(aes(x=3.709,xend=4.28,y=175,yend=175)) +
  annotate("text", x = 1.9995,y=160,label = c("Low sugar *"),family="",fontface = 3,size=4)
legend <- get_legend(T_bplot_reg)
T_bplot_reg <- T_bplot_reg + theme(legend.position="none")
T_bplot_reg

## Exploitation phase
TR_bplot_reg = ggplot(aes(y=early_explo,x=region,fill=treatment),data=subset(exp_raw_data, phase == "treward")) + theme_bw() +
  geom_boxplot(notch=F,position=position_dodge(0.88)) + ylim(0,380) +
scale_fill_manual(name="Treatment", labels = c("High sugar","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="b. Exploitation of a new food source",x="",y="Foraging activity \n (in number of workers foraging)\n") +
  scale_x_discrete(labels=c("Argentina","California","Australia", "New Zealand")) +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,0,10),"pt")) +
  annotate("text", x = c(0.906,1.093,1.28,0.9995,1.9995,2.9995,3.9995),
           y=c(rep(270,3),rep(375,4)),
           label = c("***","***","***","a","a","b","a"),fontface = 2,size=6) +
  geom_segment(aes(x=0.709,xend=1.28,y=350,yend=350)) +
  geom_segment(aes(x=1.709,xend=2.28,y=350,yend=350)) +
  geom_segment(aes(x=2.709,xend=3.28,y=350,yend=350)) +
  geom_segment(aes(x=3.709,xend=4.28,y=350,yend=350)) +
  annotate("text", x = 2.9995,y=315,label = c("OA *\nOA + Low sugar **"),family="",fontface = 3,size=4)
TR_bplot_reg

grid.arrange(T_bplot_reg,TR_bplot_reg,legend,ncol=1,nrow=3,layout_matrix=rbind(1,2,3),heights = c(2.5,2.5, 0.75))
dev.copy2pdf(file=paste(FigDir,"/fig_treatments.pdf",sep=""), 
             width=12, 
             height=8)


### FIGURE 4: Time to the top of the tree

## Individual plots per region
surv_AR <- subset(survfit_coordinates, region == "AR")
surv_CA <- subset(survfit_coordinates, region == "CA")
surv_OZ <- subset(survfit_coordinates, region == "AU")
surv_NZ <- subset(survfit_coordinates, region == "NZ")

gg_survival_AR = ggplot(aes(y=prop,x=time,fill=treatment), data=surv_AR) + theme_bw() +
  geom_point() + geom_line(data=surv_AR,aes(linetype=treatment)) +
  scale_linetype_manual(name = "Treatment",values = c(control = "solid", low_sugar = "dotted", OA = "dashed", OA_low_sugar = "dotdash"),
                        labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA")) +
  geom_ribbon(data=surv_AR,aes(ymin=lower,ymax=upper),alpha=0.5)+
  scale_fill_manual(name="95% confidence interval", labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="a. Argentina",x="",y="Proportion of colonies\nthat reached the top branch\n") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="bottom", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,2,10),"pt")) +
  annotate("text", x = rep(45,3), y=c(0.20,0.30,0.40), color= c(blue,green,red),
           label = c("Low sugar *","High sugar + OA *","Low sugar + OA ***"),fontface = 2,size=6)
legend <- get_legend(gg_survival_AR)
gg_survival_AR <- gg_survival_AR + theme(legend.position="none")
gg_survival_AR


gg_survival_CA = ggplot(aes(y=prop,x=time,fill=treatment),data=surv_CA) + theme_bw() +
  geom_point() + geom_line(data=surv_CA,aes(linetype=treatment)) +
  scale_linetype_manual(name = "Treatment",values = c(control = "solid", low_sugar = "dotted", OA = "dashed", OA_low_sugar = "dotdash"),
                        labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA")) +
  geom_ribbon(data=surv_CA,aes(ymin=lower,ymax=upper),alpha=0.5)+
  scale_fill_manual(name="Treatment", labels = c("High sugar","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="b. California",x="",y="\n") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,2,10),"pt")) +
  annotate("text", x = rep(45,2), y=c(0.35,0.20), color= c(black,blue),
           label = c("Significant interaction\ntreatment x region:","Low sugar ***"),fontface = 2,size=6)
gg_survival_CA

gg_survival_OZ = ggplot(aes(y=prop,x=time,fill=treatment), data=surv_OZ) + theme_bw() +
  geom_point() + geom_line(data=surv_OZ,aes(linetype=treatment)) +
  scale_linetype_manual(name = "Treatment",values = c(control = "solid", low_sugar = "dotted", OA = "dashed", OA_low_sugar = "dotdash"),
                        labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA")) +
  geom_ribbon(data=surv_OZ,aes(ymin=lower,ymax=upper),alpha=0.5)+
  scale_fill_manual(name="95% confidence interval", labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="c. Australia",x="\nTime after the start of the assay (in minutes)",y="Proportion of colonies\nthat reached the top branch\n") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,2,10),"pt")) +
  annotate("text", x = rep(45,1), y=c(0.30), color= c(black),
           label = c("No significant\ninteraction"),fontface = 2,size=6)
gg_survival_OZ

gg_survival_NZ = ggplot(aes(y=prop,x=time,fill=treatment), data=surv_NZ) + theme_bw() +
  geom_point() + geom_line(data=surv_NZ,aes(linetype=treatment)) +
  scale_linetype_manual(name = "Treatment",values = c(control = "solid", low_sugar = "dotted", OA = "dashed", OA_low_sugar = "dotdash"),
                        labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA")) +
  geom_ribbon(data=surv_NZ,aes(ymin=lower,ymax=upper),alpha=0.5)+
  scale_fill_manual(name="95% confidence interval", labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="d. New Zealand",x="\nTime after the start of the assay (in minutes)",y="\n") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,2,10),"pt")) +
  annotate("text", x = rep(45,1), y=c(0.30), color= c(black),
           label = c("No significant\ninteraction"),fontface = 2,size=6)
gg_survival_NZ

grid.arrange(gg_survival_AR,gg_survival_CA,gg_survival_OZ,gg_survival_NZ,legend,ncol=2,nrow=3,layout_matrix=rbind(c(1,2),c(3,4),c(5,5)),heights = c(2.5,2.5, 0.75))

dev.copy2pdf(file=paste(FigDir,"/fig_survival.pdf",sep=""), 
             width=12, 
             height=12)

### FIGURE 4: Aggression

aggression_full$bin <- ifelse(aggression_full$agg_full < 2, "Non\naggressive", "Aggressive")
aggression_full$bin <- factor(aggression_full$bin,levels=c("Non\naggressive","Aggressive",ordered=T))
aggression_full$label <- with(aggression_full,
                              ifelse(treatment == "control","Control\n(high sugar)",
                                     ifelse(treatment == "low_sugar","Low Sugar",
                                            ifelse(treatment == "OA", "High sugar\n+ OA", "Low sugar\n+ OA"))))
aggression_full$label <- factor(aggression_full$label,levels=c("Control\n(high sugar)","Low Sugar","High sugar\n+ OA", "Low sugar\n+ OA"),ordered=T)
aggression_full$region_name <- ifelse(aggression_full$region == "AR", "a. Argentina","b. Australia (p < 0.001)")

gg_agg_all <- ggplot(aggression_full, aes(x=label,fill=as.factor(agg_full))) + theme_bw() +
  geom_bar(position="fill") + ylim(0,1.04) +
  scale_fill_manual(values=c("darkblue","Dodgerblue","red4","orangered3","red"),name="",labels=c("Avoidance","Neutral","Reactive attack","Proactive Attack","Prolonged fight")) +
  facet_grid(~region_name) + labs(title="",x="",y="Proportion of interactions recorded\n") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0),strip.text = element_text(size=14), panel.grid.major.y = element_line(color = "black")) +
  theme(legend.position="bottom", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(),
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,2,10),"pt"))
gg_agg_all

dev.copy2pdf(file=paste(FigDir,"/fig_aggression_stacked.pdf",sep=""), 
             width=10, 
             height=8)

### FIGURES S2/S5: OA head titration/TYR, 5-HT titration

exp_samples$Treatment <- factor(exp_samples$Treatment,levels = c("CTRL","LS","OA","LSOA"),ordered=T)

mean_amounts <- data.frame("Treatment" = c("CTRL","LS","OA","LSOA"), "Mean" = c(17.004,11.927,87.22,226.17))

plot_OA <- ggplot(exp_samples, aes(x=Treatment,y=amount_OA,fill=Treatment)) + theme_bw() + scale_y_log10() +
  geom_boxplot(notch=F,position=position_dodge(0.88)) + geom_point(data=mean_amounts,aes(x=Treatment,y=Mean),shape=21,fill="white",color="black",size=4) +
  scale_x_discrete(labels=c("Control\n(high sugar)","Low sugar","High sugar + OA", "low sugar + OA")) +
  scale_fill_manual(name="Treatment", labels = c("Control\n(high sugar)","Low sugar","High sugar\n+ OA","Low sugar\n+ OA"),values=c(black,blue,green,red)) +
  labs(title="OA content in relation to treatment",x="",y="OA content per worker head (pg)") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="bottom", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,0,10),"pt")) +
  annotate("text", x = c(2,3,4), y=c(60,2500,2500), label = c("NS", "*", "*"),fontface = 2,size=6)
plot_OA

dev.copy2pdf(file=paste(FigDir,"/HPLC_OA.pdf",sep=""), 
             width=10, 
             height=6)

plot_SER <- ggplot(exp_samples, aes(x=Treatment,y=amount_SER,fill=Treatment)) + theme_bw() +
  geom_boxplot(notch=F,position=position_dodge(0.88)) + geom_point(data=mean_amounts,aes(x=Treatment,y=Mean),shape=21,fill="white",color="black",size=4) +
  scale_x_discrete(labels=c("Control\n(high sugar)","Low sugar","High sugar + OA", "low sugar + OA")) +
  scale_fill_manual(name="Treatment", labels = c("Control\n(high sugar)","Low sugar","High sugar\n+ OA","Low sugar\n+ OA"),values=c(black,blue,green,red)) +
  labs(title="a. 5-HT content in relation to treatment",x="",y="5-HT content per worker head (pg)") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="none", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,0,10),"pt"))# +
plot_SER

plot_TYR <- ggplot(exp_samples, aes(x=Treatment,y=amount_TYR,fill=Treatment)) + theme_bw() +scale_y_log10() +
  geom_boxplot(notch=F,position=position_dodge(0.88)) + geom_point(data=mean_amounts,aes(x=Treatment,y=Mean),shape=21,fill="white",color="black",size=4) +
  scale_x_discrete(labels=c("Control\n(high sugar)","Low sugar","High sugar + OA", "low sugar + OA")) +
  scale_fill_manual(name="Treatment", labels = c("Control\n(high sugar)","Low sugar","High sugar\n+ OA","Low sugar\n+ OA"),values=c(black,blue,green,red)) +
  labs(title="b. TYR content in relation to treatment",x="",y="TYR content per worker head (pg)") +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="bottom", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,0,10),"pt")) #+
legend <- get_legend(plot_TYR)
plot_TYR <- plot_TYR + theme(legend.position="none")
plot_TYR

grid.arrange(plot_SER,plot_TYR,legend,ncol=1,nrow=3,layout_matrix=rbind(1,2,3),heights = c(2.5,2.5,0.5))

dev.copy2pdf(file=paste(FigDir,"/HPLC_SERTYR.pdf",sep=""), 
width=10, 
height=12)
dev.off()

### FIGURE S3: Foraging profiles

par(mfrow = c(4, 8)) # 2-by-2 grid of plots
par(oma = c(4, 4, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = rep(2.75, 4)) # make the plots be closer together

exp_raw_data_names <- exp_raw_data
exp_raw_data_names$treatment <- as.character(exp_raw_data_names$treatment); exp_raw_data_names$phase <- as.character(exp_raw_data_names$phase)
exp_raw_data_names$treatment[exp_raw_data_names$treatment == "control"] <- as.character("Control")
exp_raw_data_names$treatment[exp_raw_data_names$treatment == "low_sugar"] <- as.character("Low sugar")
exp_raw_data_names$treatment[exp_raw_data_names$treatment == "OA_low_sugar"] <- as.character("Low sugar + OA")
exp_raw_data_names$phase[exp_raw_data_names$phase == "treatment"] <- as.character("Exploration")
exp_raw_data_names$phase[exp_raw_data_names$phase == "treward"] <- as.character("Exploitation")
exp_raw_data_names$treatment <- as.factor(exp_raw_data_names$treatment); exp_raw_data_names$phase <- as.factor(exp_raw_data_names$phase)
phase_list = c("Exploration", "Exploitation")
treatment_list = c("Control", "Low sugar", "OA", "Low sugar + OA")
region_list = c("AR","CA","AU","NZ")
x = c(0, 10, 20, 30, 40, 50, 60)
slopes = NULL
for (country in region_list) {
  data_region = subset(exp_raw_data_names, region == country)
  for (exp_phase in phase_list) {
    data_phase = subset(data_region, phase == exp_phase)
    for (condition in treatment_list) {
      
      data = subset(data_phase, treatment == condition)
      y = c(mean(data[, 11]), mean(data[, 12]), mean(data[, 13]), mean(data[, 14]), mean(data[, 15]), 
            mean(data[, 16]), mean(data[, 17]))
      sd = c(sd(data[, 11]), sd(data[, 12]), sd(data[, 13]), sd(data[, 14]), sd(data[, 15]), sd(data[, 
                                                                                                     16]), sd(data[, 17]))
      plot(x, y, xlim = c(0, 65), ylim = c(-10, 120), type = "n", xlab= "",ylab="", bty = "l", main = c(country,exp_phase, condition))
      with(data, expr = errbar(x, y, y + sd, y - sd, add = T, pch = 16, cap = 0.05, type = "b", bty = "l"))
      
      sum = summary(lm(y[1:4] ~ x[1:4]))
      line_slope = data.frame(slope = sum$coefficients[2, 1], se = sum$coefficients[2, 2], phase = exp_phase, 
                              treatment = condition, region = country)
      slopes = rbind(slopes, line_slope)
    }
  }
}

# print the overall labels
mtext('Time (minutes)', side = 1, outer = TRUE, line = 2,cex=1.5)
mtext('Foraging activity (number of workers foraging)', side = 2, outer = TRUE, line = 2,cex=1.5)

dev.copy2pdf(file=paste(FigDir,"/foraging_profiles.pdf",sep=""), 
             width=14, 
             height=8)


### FIGURE S4: Effect of treatmnets along the introduction pathway on sugar discovery

bplot_sug = ggplot(aes(y=sug_10,x=region,fill=treatment),data=subset(exp_raw_data, phase == "treward")) + geom_boxplot(notch=F,position=position_dodge(0.88)) + theme_bw() +
  scale_fill_manual(name="Treatment", labels = c("Control\n(high sugar)","Low sugar","High sugar + OA","Low sugar + OA"),values=c(black,blue,green,red)) +
  labs(title="Foraging performance associated with resource discovery",x="",y="Number of workers feeding\non a new food souce after 10 minutes\n") +
  scale_x_discrete(labels=c("Argentina","California","Australia", "New Zealand")) +
  theme(plot.title = element_text(face="bold", size=18, hjust=0)) +
  theme(legend.position="bottom", legend.title=element_text(face="bold"), legend.text=element_text(size=14), 
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.text.x = element_text(size = 14, colour = "black"), 
        axis.title.y=element_text(size = 16,vjust=1, colour = "black"),
        axis.title.x=element_text(size = 16, colour = "black"),
        plot.margin=unit(c(2,10,0,10),"pt")) +
  annotate("text", x = c(0.906,1.093,1.28,0.9995,1.9995,2.9995,3.9995),
           y=c(rep(80,3),rep(100,4)),
           label = c("***","***","***","a","a","a","b"),fontface = 2,size=6) +
  geom_segment(aes(x=0.709,xend=1.28,y=95,yend=95)) +
  geom_segment(aes(x=1.709,xend=2.28,y=95,yend=95)) +
  geom_segment(aes(x=2.709,xend=3.28,y=95,yend=95)) +
  geom_segment(aes(x=3.709,xend=4.28,y=95,yend=95))
bplot_sug
dev.copy2pdf(file=paste(FigDir,"/fig_sugar.pdf",sep=""), 
             width=10, 
             height=6)
