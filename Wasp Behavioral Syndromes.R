###### Code to analyze wasp behavioral syndromes
###### Written by Emily Laub (eclaub@umich.edu) in 2023

################## Packages for analysis

library(dplyr)
library(ggplot2)
library(rptR)
library(MCMCglmm)

################# load data

wasp_personality <- file.choose("Behavioral Syndromes for submission 20231106.csv")
personality <- read.csv(wasp_personality)

############################################ Graphing template
# Starts with theme_grey and then modify some parts
theme_mine <- function(base_size = 18, base_family = "Helvetica") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.text.x = element_text(size=20),
      axis.text.y = element_text(size=20,hjust=1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size=25),
      axis.title.y= element_text(size=25,angle=90),
      panel.background = element_blank(), 
      panel.border =element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line.x = element_line(color="black", size = 1),
      axis.line.y = element_line(color="black", size = 1)
    )
}
################################################################
################## Repeatability of personality traits

### Activity, only ID
activity1 <- rpt(seconds_moving ~ (1|Wasp.ID), grname = c("Wasp.ID"), data = personality,
                 datatype = "Gaussian", nboot = 1000, ratio = TRUE)
summary(activity1)

### Activity, Fixed Effects
activity2 <- rpt(seconds_moving ~ Trial + weight + (1|Wasp.ID), grname = c("Wasp.ID", "Fixed", "Residual"), data = explore,
                 datatype = "Gaussian", nboot = 1000, ratio = TRUE)
summary(activity2)

activity3 <- rpt(seconds_moving ~ Trial + weight + (1|Wasp.ID), grname = c("Wasp.ID", "Fixed", "Residual"), data = explore,
                 datatype = "Gaussian", nboot = 1000, ratio = FALSE)
summary(activity3)

### Exploration, only ID
exploration1 <- rpt(chambers_entered ~ (1|Wasp.ID), grname = c("Wasp.ID"), data = personality,
                 datatype = "Poisson", nboot = 1000, ratio = TRUE)
summary(exploration1)

### Exploration, Fixed Effects
exploration2 <- rpt(chambers_entered ~ Trial + weight + (1|Wasp.ID), grname = c("Wasp.ID", "Fixed", "Residual"), data = personality,
                 datatype = "Poisson", nboot = 1000, ratio = TRUE)
summary(exploration2)

exploration3 <- rpt(chambers_entered ~ Trial + weight + (1|Wasp.ID), grname = c("Wasp.ID", "Fixed", "Residual"), data = personality,
                 datatype = "Poisson", nboot = 1000, ratio = FALSE)
summary(exploration3)

### Affiliation, only ID
aff1<- rpt(bodily_contact_time ~ (1|Wasp.ID), 
           grname =  c("Wasp.ID"), data = personality,
           datatype = "Gaussian", nboot = 1000, ratio = TRUE)
summary(aff1)

### Affiliation, Fixed Effects
aff2 <- rpt(bodily_contact_time ~ Dummy_Color + Trial + weight  + (1|Wasp.ID), 
           grname =  c("Wasp.ID", "Fixed", "Residual"), data = personality,
           datatype = "Gaussian", nboot = 1000, ratio = TRUE)
summary(aff2)

aff3 <- rpt(bodily_contact_time ~ Dummy_Color + Trial + weight  + (1|Wasp.ID), 
           grname =  c("Wasp.ID", "Fixed", "Residual"), data = personality,
           datatype = "Gaussian", nboot = 1000, ratio = FALSE)
summary(aff3)

### Aggression, only ID
agg1 <- rpt(Log_aggressionression ~ (1|Wasp.ID), grname = "Wasp.ID", data = personality, datatype = "Gaussian",
            nboot = 1000)
summary(agg1)

### Aggression, Fixed effects
agg2 <- rpt(Log_aggressionression ~ Dummy_Color + Trial + weight + (1|Wasp.ID), 
            grname = c("Wasp.ID", "Fixed", "Residual"), data = personality,
            nboot= 1000, datatype = "Gaussian", ratio = TRUE)
summary(agg2)

agg3 <- rpt(Log_aggressionression ~ Dummy_Color + Trial + weight + (1|Wasp.ID), 
            grname = c("Wasp.ID", "Fixed", "Residual"), data = personality,
            nboot= 1000, datatype = "Gaussian", ratio = FALSE)
summary(agg3)

### Neutral behavior, only ID
neutral1 <- rpt(anntenation~ (1|Wasp.ID), grname = "Wasp.ID", data = personality, datatype = "Poisson",
            nboot = 1000)
summary(neutral1)

### Neutral behavior, Fixed effects
neutral2 <- rpt(anntenation~ Dummy_Color + Trial + weight + (1|Wasp.ID), 
                grname = c("Wasp.ID", "Fixed", "Residual"), data = personality,
                nboot= 1000, datatype = "Poisson", ratio = TRUE)
summary(neutral2)

neutral3 <- rpt(anntenation~ Dummy_Color + Trial + weight + (1|Wasp.ID), 
                grname = c("Wasp.ID", "Fixed", "Residual"), data = personality,
                nboot= 1000, datatype = "Poisson", ratio = FALSE)
summary(neutral3)


#########################################################################################################
################################################## Correlations between behaviors

############################# Spearman's rank correlations
####### Get averages across trials
personality_summary =
  personality %>% 
  group_by(Wasp.ID) %>%
  summarise(Affiliation = mean(bodily_contact_time),
            Aggression = mean(Log_aggressionression), Neutral = mean(anntenation), Activity = mean(seconds_moving), 
            Exploration = mean(chambers_entered)) %>%
  as.data.frame()

####### Spearman's rank correlations
cor.test(summary_all_pers$Affiliation, summary_all_pers$Activity, method = c("spearman"))
cor.test(summary_all_pers$Affiliation, summary_all_pers$Aggression, method = c("spearman"))
cor.test(summary_all_pers$Affiliation, summary_all_pers$Exploration, method = c("spearman"))
cor.test(summary_all_pers$Affiliation, summary_all_pers$Neutral, method = c("spearman"))

cor.test(summary_all_pers$Aggression, summary_all_pers$Activity, method = c("spearman"))
cor.test(summary_all_pers$Aggression, summary_all_pers$Exploration, method = c("spearman"))
cor.test(summary_all_pers$Aggression, summary_all_pers$Neutral, method = c("spearman"))

cor.test(summary_all_pers$Activity, summary_all_pers$Exploration, method = c("spearman"))
cor.test(summary_all_pers$Activity, summary_all_pers$Neutral, method = c("spearman"))

cor.test(summary_all_pers$Exploration, summary_all_pers$Neutral, method = c("spearman"))

################ Graphing significant correlations

ggplot(summary_all_pers, aes(Activity, Exploration)) +
  geom_point() +
  theme_mine() +
  xlab("Activity") + ylab("Exploration") +
  geom_smooth(method = "lm") 

ggplot(summary_all_pers, aes(Activity, Aggression)) +
  geom_point() +
  theme_mine() +
  xlab("Activity") + ylab("Aggression") +
  geom_smooth(method = "lm") 

ggplot(summary_all_pers, aes(Activity, Affiliation)) +
  geom_point() +
  theme_mine() +
  xlab("Activity") + ylab("Affiliation") +
  geom_smooth(method = "lm") 

############ MCMC 

personality$chambers_entered <- as.numeric(personality$chambers_entered)
personality$Trial <- as.numeric(personality$Trial)
prior_A_E_1px = list(R = list(V = diag(2), nu = 0.002),
                     G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                        alpha.V = diag(25^2,2,2))))

mcmc_A_E_us <- MCMCglmm(cbind(scale(seconds_moving), (chambers_entered)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                          trait:scale(weight),
                        random =~ us(trait):Wasp.ID,
                        rcov =~ us(trait):units,
                        family = c("gaussian","poisson"), prior = prior_A_E_1px, nitt=420000,
                        burnin=20000,
                        thin=100,
                        verbose = TRUE,
                        data = as.data.frame(personality))

summary(mcmc_A_E_us)
mcmc_cor_AE <- mcmc_A_E_us$VCV[,"traitseconds_moving:traitchambers_entered.Wasp.ID"]/ (sqrt(mcmc_A_E_us$VCV[,"traitseconds_moving:traitseconds_moving.Wasp.ID"])*
                                                                                             sqrt(mcmc_A_E_us$VCV[,"traitchambers_entered:traitchambers_entered.Wasp.ID"]))

mean(mcmc_cor_AE)

HPDinterval(mcmc_cor_AE)



##############################################
############ Activity and Affiliation
head(personality)
str(personality)
personality$bodily_contact_time <- as.numeric(personality$bodily_contact_time)
prior_A_A_1px = list(R = list(V = diag(2), nu = 0.002),
                     G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                        alpha.V = diag(25^2,2,2))))

mcmc_A_A_us <- MCMCglmm(cbind(scale(seconds_moving), scale(bodily_contact_time)) ~ trait-1 + at.level(trait,2):(as.numeric(as.factor(Dummy_Color))) + trait:scale(Trial, scale = FALSE) +
                          trait:scale(weight),
                        random =~ us(trait):Wasp.ID,
                        rcov =~ us(trait):units,
                        family = c("gaussian","gaussian"), prior = prior_A_A_1px, nitt=420000,
                        burnin=20000,
                        thin=100,
                        verbose = TRUE,
                        data = as.data.frame(personality))

summary(mcmc_A_A_us)
mcmc_cor_AA <- mcmc_A_A_us$VCV[,"traitseconds_moving:traitbodily_contact_time.Wasp.ID"]/ (sqrt(mcmc_A_A_us$VCV[,"traitseconds_moving:traitseconds_moving.Wasp.ID"])*
                                                                                      sqrt(mcmc_A_A_us$VCV[,"traitbodily_contact_time:traitbodily_contact_time.Wasp.ID"]))

mean(mcmc_cor_AA)

HPDinterval(mcmc_cor_AA)


################ Aggressive and Active
prior_Ac_Agg_1px = list(R = list(V = diag(2), nu = 0.002),
                        G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                           alpha.V = diag(25^2,2,2))))
mcmc_Ac_Agg_us <- MCMCglmm(cbind(scale(seconds_moving), scale(Log_aggression)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                             trait:scale(weight) + at.level(trait,2):(as.numeric(as.factor(Dummy_Color))),
                           random =~ us(trait):Wasp.ID,
                           rcov =~ us(trait):units,
                           family = c("gaussian","gaussian"), prior = prior_Ac_Agg_1px, nitt=420000,
                           burnin=20000,
                           thin=100,
                           verbose = TRUE,
                           data = as.data.frame(personality))
summary(mcmc_Ac_Agg_us)
mcmc_cor_AcAgg <- mcmc_Ac_Agg_us$VCV[,"traitseconds_moving:traitLog_aggression.Wasp.ID"]/ (sqrt(mcmc_Ac_Agg_us$VCV[,"traitseconds_moving:traitseconds_moving.Wasp.ID"])*
                                                                                sqrt(mcmc_Ac_Agg_us$VCV[,"traitLog_aggression:traitLog_aggression.Wasp.ID"]))


mean(mcmc_cor_AcAgg)
HPDinterval(mcmc_cor_AcAgg)


########################### Aggression and Exploratory.
prior_Ex_Agg_1px = list(R = list(V = diag(2), nu = 0.002),
                        G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                           alpha.V = diag(25^2,2,2))))


mcmc_Ex_Agg_us <- MCMCglmm(cbind((chambers_entered), scale(Log_aggression)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                             trait:scale(weight) + at.level(trait,2):(as.numeric(as.factor(Dummy_Color))),
                           random =~ us(trait):Wasp.ID,
                           rcov =~ us(trait):units,
                           family = c("poisson","gaussian"), prior = prior_Ex_Agg_1px, nitt=420000,
                           burnin=20000,
                           thin=100,
                           verbose = TRUE,
                           data = as.data.frame(personality))
summary(mcmc_Ex_Agg_us)
mcmc_cor_ExAgg <- mcmc_Ex_Agg_us$VCV[,"traitchambers_entered:traitLog_aggression.Wasp.ID"]/ (sqrt(mcmc_Ex_Agg_us$VCV[,"traitchambers_entered:traitchambers_entered.Wasp.ID"])*
                                                                                                  sqrt(mcmc_Ex_Agg_us$VCV[,"traitLog_aggression:traitLog_aggression.Wasp.ID"]))

mean(mcmc_cor_ExAgg)
HPDinterval(mcmc_cor_ExAgg)


############################# Activity and Neutral
prior_Ac_Neu_1px = list(R = list(V = diag(2), nu = 0.002),
                        G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                           alpha.V = diag(25^2,2,2))))

mcmc_Ac_Neu_us <- MCMCglmm(cbind(scale(seconds_moving), (anntenation)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                             trait:scale(weight) + at.level(trait,2):(as.numeric(as.factor(Dummy_Color))),
                           random =~ us(trait):Wasp.ID,
                           rcov =~ us(trait):units,
                           family = c("gaussian","poisson"), prior = prior_Ac_Neu_1px, nitt=420000,
                           burnin=20000,
                           thin=100,
                           verbose = TRUE,
                           data = as.data.frame(personality))
summary(mcmc_Ac_Neu_us)
mcmc_cor_AcNeu <- mcmc_Ac_Neu_us$VCV[,"traitseconds_moving:traitanntenation.Wasp.ID"]/ (sqrt(mcmc_Ac_Neu_us$VCV[,"traitseconds_moving:traitseconds_moving.Wasp.ID"])*
                                                                                      sqrt(mcmc_Ac_Neu_us$VCV[,"traitanntenation:traitanntenation.Wasp.ID"]))

mean(mcmc_cor_AcNeu)
HPDinterval(mcmc_cor_AcNeu)

######################## Exploration and Neutral

prior_Ex_Neu_1px = list(R = list(V = diag(2), nu = 0.002),
                        G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                           alpha.V = diag(25^2,2,2))))

mcmc_Ex_Neu_us <- MCMCglmm(cbind((chambers_entered), (anntenation)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                             trait:scale(weight)+ at.level(trait,2):(as.numeric(as.factor(Dummy_Color))),
                           random =~ us(trait):Wasp.ID,
                           rcov =~ us(trait):units,
                           family = c("poisson","poisson"), prior = prior_Ex_Neu_1px, nitt=420000,
                           burnin=20000,
                           thin=100,
                           verbose = TRUE,
                           data = as.data.frame(personality))
summary(mcmc_Ex_Neu_us)
mcmc_cor_ExNeu <- mcmc_Ex_Neu_us$VCV[,"traitchambers_entered:traitanntenation.Wasp.ID"]/ (sqrt(mcmc_Ex_Neu_us$VCV[,"traitchambers_entered:traitchambers_entered.Wasp.ID"])*
                                                                                                        sqrt(mcmc_Ex_Neu_us$VCV[,"traitanntenation:traitanntenation.Wasp.ID"]))

mean(mcmc_cor_ExNeu)
HPDinterval(mcmc_cor_ExNeu)

############################## Aggression and Neutral
prior_Agg_Neu_1px = list(R = list(V = diag(2), nu = 0.002),
                         G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                            alpha.V = diag(25^2,2,2))))

mcmc_Agg_Neu_us <- MCMCglmm(cbind(scale(Log_aggression), (anntenation)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                              trait:scale(weight) + at.level(trait,2):(as.numeric(Dummy_Color)) + at.level(trait,1):(as.numeric(as.factor(Dummy_Color))),
                            random =~ us(trait):Wasp.ID,
                            rcov =~ us(trait):units,
                            family = c("gaussian","poisson"), prior = prior_Agg_Neu_1px, nitt=420000,
                            burnin=20000,
                            thin=100,
                            verbose = TRUE,
                            data = as.data.frame(personality))
summary(mcmc_Agg_Neu_us)
mcmc_cor_AggNeu <- mcmc_Agg_Neu_us$VCV[,"traitLog_aggression:traitanntenation.Wasp.ID"]/ (sqrt(mcmc_Agg_Neu_us$VCV[,"traitLog_aggression:traitLog_aggression.Wasp.ID"])*
                                                                                       sqrt(mcmc_Agg_Neu_us$VCV[,"traitanntenation:traitanntenation.Wasp.ID"]))

mean(mcmc_cor_AggNeu)
HPDinterval(mcmc_cor_AggNeu)

########################## Affiliation and Neutral

prior_Aff_Neu_1px = list(R = list(V = diag(2), nu = 0.002),
                         G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                            alpha.V = diag(25^2,2,2))))

mcmc_Aff_Neu_us <- MCMCglmm(cbind(scale(bodily_contact_time), (anntenation)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                              trait:scale(weight) + at.level(trait,2):(as.numeric(Dummy_Color)) + at.level(trait,1):(as.numeric(as.factor(Dummy_Color))),
                            random =~ us(trait):Wasp.ID,
                            rcov =~ us(trait):units,
                            family = c("gaussian","poisson"), prior = prior_Aff_Neu_1px, nitt=420000,
                            burnin=20000,
                            thin=100,
                            verbose = TRUE,
                            data = as.data.frame(personality))
summary(mcmc_Aff_Neu_us)
mcmc_cor_AffNeu <- mcmc_Aff_Neu_us$VCV[,"traitbodily_contact_time:traitanntenation.Wasp.ID"]/ (sqrt(mcmc_Aff_Neu_us$VCV[,"traitbodily_contact_time:traitbodily_contact_time.Wasp.ID"])*
                                                                                                   sqrt(mcmc_Aff_Neu_us$VCV[,"traitanntenation:traitanntenation.Wasp.ID"]))

mean(mcmc_cor_AffNeu)
HPDinterval(mcmc_cor_AffNeu)


######################## Affiliation and Aggression

prior_Aff_Agg_1px = list(R = list(V = diag(2), nu = 0.002),
                         G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                            alpha.V = diag(25^2,2,2))))

mcmc_Aff_Agg_us <- MCMCglmm(cbind(scale(bodily_contact_time), scale(Log_aggression)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                              trait:scale(weight) + at.level(trait,2):(as.numeric(as.factor(Dummy_Color))) + at.level(trait,1):(as.numeric(as.factor(Dummy_Color))),
                            random =~ us(trait):Wasp.ID,
                            rcov =~ us(trait):units,
                            family = c("gaussian","gaussian"), prior = prior_Aff_Agg_1px, nitt=420000,
                            burnin=20000,
                            thin=100,
                            verbose = TRUE,
                            data = as.data.frame(personality))
summary(mcmc_Aff_Agg_us)
mcmc_cor_AffAgg <- mcmc_Aff_Agg_us$VCV[,"traitbodily_contact_time:traitLog_aggression.Wasp.ID"]/ (sqrt(mcmc_Aff_Agg_us$VCV[,"traitbodily_contact_time:traitbodily_contact_time.Wasp.ID"])*
                                                                                             sqrt(mcmc_Aff_Agg_us$VCV[,"traitLog_aggression:traitLog_aggression.Wasp.ID"]))

mean(mcmc_cor_AffAgg)
HPDinterval(mcmc_cor_AffAgg)



mcmc_prop_B <- mcmc_Aff_Agg_us$VCV[,"traitLog_aggression:traitLog_aggression.Wasp.ID"]/(
  mcmc_Aff_Agg_us$VCV[,"traitLog_aggression:traitLog_aggression.Wasp.ID"] +
    mcmc_Aff_Agg_us$VCV[,"traitLog_aggression:traitLog_aggression.units"]
)
mean(mcmc_prop_B)
HPDinterval(mcmc_prop_B)

########################## Affiliation and Exploration 

prior_Aff_Ex_1px = list(R = list(V = diag(2), nu = 0.002),
                        G = list(G1 = list(V = diag(2), nu = 2, alpha.mu = rep(0,2),
                                           alpha.V = diag(25^2,2,2))))

mcmc_Aff_Ex_us <- MCMCglmm(cbind(scale(bodily_contact_time), (chambers_entered)) ~ trait-1 + trait:scale(Trial, scale = FALSE) +
                             trait:scale(weight) + at.level(trait,1):scale(as.numeric(as.factor(Dummy_Color))),
                           random =~ us(trait):Wasp.ID,
                           rcov =~ us(trait):units,
                           family = c("gaussian","poisson"), prior = prior_Aff_Ex_1px, nitt=420000,
                           burnin=20000,
                           thin=100,
                           verbose = TRUE,
                           data = as.data.frame(personality))
summary(mcmc_Aff_Ex_us)
mcmc_cor_AffEx <- mcmc_Aff_Ex_us$VCV[,"traitbodily_contact_time:traitchambers_entered.Wasp.ID"]/ (sqrt(mcmc_Aff_Ex_us$VCV[,"traitbodily_contact_time:traitbodily_contact_time.Wasp.ID"])*
                                                                                                              sqrt(mcmc_Aff_Ex_us$VCV[,"traitchambers_entered:traitchambers_entered.Wasp.ID"]))

mean(mcmc_cor_AffEx)
HPDinterval(mcmc_cor_AffEx)



################# Graph relationship below

df_syndromes <- data_frame(Traits = c("Activity, Exploration", 
                                      "Activity, Affiliation",
                                      "Activity, Aggression", 
                                      "Activity, Antennation", 
                                      "Exploration, Affiliation", 
                                      "Exploration, Aggression", 
                                      "Exploration, Antennation", 
                                      "Affiliation, Aggression", 
                                      "Affiliation, Antennation", 
                                      "Aggression, Antennation"), 
                           
                           Estimate = c(mean(mcmc_cor_AE_Poisson),
                                        mean(mcmc_cor_AA),
                                        mean(mcmc_cor_AcAgg),
                                        mean(mcmc_cor_AcNeu),
                                        mean(mcmc_cor_AffEx),
                                        mean(mcmc_cor_ExAgg),
                                        mean(mcmc_cor_ExNeu),
                                        mean(mcmc_cor_AffAgg),
                                        mean(mcmc_cor_AffNeu),
                                        mean(mcmc_cor_AggNeu)),
                           
                           Lower = c(HPDinterval(mcmc_cor_AE_Poisson)[,"lower"],
                                     HPDinterval(mcmc_cor_AA)[,"lower"],
                                     HPDinterval(mcmc_cor_AcAgg)[,"lower"],
                                     HPDinterval(mcmc_cor_AcNeu)[,"lower"],
                                     HPDinterval(mcmc_cor_AffEx)[,"lower"],
                                     HPDinterval(mcmc_cor_ExAgg)[,"lower"],
                                     HPDinterval(mcmc_cor_ExNeu)[,"lower"],
                                     HPDinterval(mcmc_cor_AffAgg)[,"lower"],
                                     HPDinterval(mcmc_cor_AffNeu)[,"lower"],
                                     HPDinterval(mcmc_cor_AggNeu)[,"lower"]
                           ), 
                           
                           Upper = c(HPDinterval(mcmc_cor_AE_Poisson)[,"upper"],
                                     HPDinterval(mcmc_cor_AA)[,"upper"],
                                     HPDinterval(mcmc_cor_AcAgg)[,"upper"],
                                     HPDinterval(mcmc_cor_AcNeu)[,"upper"],
                                     HPDinterval(mcmc_cor_AffEx)[,"upper"],
                                     HPDinterval(mcmc_cor_ExAgg)[,"upper"],
                                     HPDinterval(mcmc_cor_ExNeu)[,"upper"],
                                     HPDinterval(mcmc_cor_AffAgg)[,"upper"],
                                     HPDinterval(mcmc_cor_AffNeu)[,"upper"],
                                     HPDinterval(mcmc_cor_AggNeu)[,"upper"]))


ggplot(df_syndromes, aes(x = Traits, y = Estimate)) + 
  geom_pointrange(aes(ymin = Lower, ymax = Upper)) + 
  geom_hline(yintercept = 0, linetype = "dotted",alpha = 0.3) + 
  scale_x_discrete(limits = c("Activity, Exploration", 
                              "Activity, Affiliation",
                              "Activity, Aggression", 
                              "Activity, Antennation", 
                              "Exploration, Affiliation", 
                              "Exploration, Aggression", 
                              "Exploration, Antennation", 
                              "Affiliation, Aggression", 
                              "Affiliation, Antennation", 
                              "Aggression, Antennation")) + 
  labs(x = "Behavior combinations", y = "Correlation (Estimate +/- 95% CIs)") + 
  ylim(-1,1) +
  coord_flip() + theme_classic(base_size = 14)




