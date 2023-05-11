# Survival analysis

##########
# See:
# https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html
# https://bioconnector.github.io/workshops/r-survival.html


# install.packages(c("survival", "survminer", "ggsurvfit", "gtsummary", "tidycmprsk"))
# devtools::install_github("zabore/condsurv")
# devtools::install_github("zabore/ezfun")

library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(condsurv)
library(survminer)
library(ggpubr)
library(viridis)

# file.choose()
setwd("~/Desktop/Vanessa survival/")

x <- "masc.survival.csv"
surv <- read.csv(x)
head(surv)
surv <- surv[,-8]
surv$individual_id <- with(surv, paste(sgRNA,oviposition_date,batch,individual_id, sep = "_"))
head(surv)

# How many individuals in the dataset?: 856
length(unique(surv$individual_id))

# How many individuals by stage and treatment?
surv %>%
  group_by(stage, sgRNA) %>%
  summarise(n=n()) %>%
  pivot_wider(names_from = "stage", values_from = "n")

# Chi-squared test for goodness of fit
# https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/
# for differences in survival to larval stage

contingency <- data.frame(
  "none" = c(306, 66), "masc284" = c(275, 19),
  row.names = c("embryos", "larvae")
)
colnames(contingency) <- c("none", "masc284")
mosaicplot(contingency, color = TRUE)
chisq.test(contingency)
# X-squared = 17.765, df = 1, p-value = 2.499e-05

contingency <- data.frame(
  "none" = c(306, 66), "masc524" = c(275, 17),
  row.names = c("embryos", "larvae")
)
colnames(contingency) <- c("none", "masc524")
mosaicplot(contingency, color = TRUE)
chisq.test(contingency)
# X-squared = 20.176, df = 1, p-value = 7.064e-06

contingency <- data.frame(
  "none" = c(306, 50), "masc284" = c(275, 15),
  row.names = c("embryos", "pupae")
)
colnames(contingency) <- c("none", "masc284")
mosaicplot(contingency, color = TRUE)
chisq.test(contingency)
# X-squared = 12.939, df = 1, p-value = 0.0003218

contingency <- data.frame(
  "none" = c(306, 50), "masc524" = c(275, 14),
  row.names = c("embryos", "pupae")
)
colnames(contingency) <- c("none", "masc524")
mosaicplot(contingency, color = TRUE)
chisq.test(contingency)
# X-squared = 14.096, df = 1, p-value = 0.0001738

contingency <- data.frame(
  "none" = c(306, 46), "masc284" = c(275, 12),
  row.names = c("embryos", "adult")
)
colnames(contingency) <- c("none", "masc284")
mosaicplot(contingency, color = TRUE)
chisq.test(contingency)
# X-squared = 14.072, df = 1, p-value = 0.0001759

contingency <- data.frame(
  "none" = c(306, 46), "masc524" = c(275, 11),
  row.names = c("embryos", "adult")
)
colnames(contingency) <- c("none", "masc524")
mosaicplot(contingency, color = TRUE)
chisq.test(contingency)
# X-squared = 15.38, df = 1, p-value = 8.791e-05


# Kaplan-Meier curves

# Hatching probability
eggs <- surv %>%
  filter(stage=="embryo") %>%
  mutate(sgRNA = sgRNA)

i <- which(eggs$status==2)
eggs$status[i] <- 0 # convert "molted" subjects to censored

# log-rank p-value
survdiff(Surv(days, status) ~ sgRNA, data = eggs)

# log-rank p-value, just using the high-quality control individuals
eggs2 <- eggs %>%
  filter(!(sgRNA=="none" & batch==1))
(x <- survdiff(Surv(days, status) ~ sgRNA, data = eggs2))
egg.pvalue <- paste0("p = ",signif(x$pvalue,2))

plot.hatching <- survfit(Surv(days, status) ~ sgRNA, data = eggs) %>%
  ggsurvplot(
    conf.int=TRUE, pval=egg.pvalue, # risk.table="abs_pct",
    title="Embryonic survival to hatching",
    legend.title="sgRNA",
    palette = "jco"
  ) +
  labs(x="days")

plot.hatching


# Pupation probability
larvae <- surv %>%
  filter(stage=="larva")

i <- which(larvae$status==2)
larvae$status[i] <- 0 # convert "molted" subjects to censored

# log-rank p-value
survdiff(Surv(days, status) ~ sgRNA, data = larvae)

# log-rank p-value, just using the high-quality control individuals
larvae2 <- larvae %>%
  filter(!(sgRNA=="none" & batch==1))
(x <- survdiff(Surv(days, status) ~ sgRNA, data = larvae2))
larvae.pvalue <- paste0("p = ",signif(x$pvalue,2))

plot.pupation <- survfit(Surv(days, status) ~ sgRNA, data = larvae) %>%
  ggsurvplot(
    conf.int=TRUE, pval=larvae.pvalue, # risk.table="abs_pct",
    title="Larval survival to pupation",
    legend.title="sgRNA",
    palette = "jco"
  ) +
  labs(x="days")

plot.pupation


# Eclosure probability
pupae <- surv %>%
  filter(stage=="pupa")

i <- which(pupae$status==2)
pupae$status[i] <- 0 # convert "molted" subjects to censored

# log-rank p-value
survdiff(Surv(days, status) ~ sgRNA, data = pupae)

# log-rank p-value, just using the high-quality control individuals
pupae2 <- pupae %>%
  filter(!(sgRNA=="none" & batch==1))
(x <- survdiff(Surv(days, status) ~ sgRNA, data = pupae2))
pupae.pvalue <- paste0("p = ",signif(x$pvalue,2))

plot.eclosure <- survfit(Surv(days, status) ~ sgRNA, data = pupae) %>%
  ggsurvplot(
    conf.int=TRUE, pval=pupae.pvalue, # risk.table="abs_pct",
    title="Pupal survival to eclosure",
    legend.title="sgRNA",
    palette = "jco"
  ) +
  labs(x="days")

plot.eclosure


# Alternative visualizations

surv.ordered <- surv %>%
  pivot_wider(id_cols = "individual_id", names_from = "stage", values_from = "days", values_fill = 0) %>%
  mutate(
    sgRNA = str_split_fixed(individual_id,"_",2)[,1],
    total = embryo + larva + pupa + adult
  ) %>%
  group_by(sgRNA) %>%
  arrange(total, embryo, larva, pupa, adult) %>%
  mutate(index = row_number()) %>%
  select(-individual_id, -total) %>%
  pivot_longer(cols = 1:4, names_to = "stage", values_to = "days") %>%
  filter(days != 0) %>%
  mutate(sgRNA = sub("none"," none",sgRNA)) %>%
  mutate(stage = sub("adult","   adult",stage)) %>%
  mutate(stage = sub("pupa","  pupa",stage)) %>%
  mutate(stage = sub("larva"," larva",stage))

surv.ordered %>%
  ggplot(aes(x=index, y=days, fill=as.factor(stage))) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y =element_blank(),
    strip.background = element_blank()
  ) +
  geom_hline(yintercept = 1:7*5, color = "gray95") +
  geom_col(width = 1) +
  scale_fill_viridis(
    name = "stage",
    discrete = TRUE, begin = 0.2, end = 0.8) +
  guides(fill = guide_legend(reverse = TRUE)) +
  facet_wrap(~sgRNA, scales="free_y") +
  ylab("days after egg-laying") +
  coord_flip()

ggsave("survival.by.individual.stage.png", width = 6, height = 6, scale = 1)
