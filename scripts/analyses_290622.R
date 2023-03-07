library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(lme4)

#analysis of bacterial densities

phpl <- read.csv("github/time_series_week_12.csv", header = T)

#calculate densities
phpl <- mutate(phpl,
               dens = count * dil * dil_fact * froz,
               dens_pls = plas_carry * dil * dil_fact * froz,
               logden = log10(dens),
               logdpl = log10(dens_pls),
               plas_prop = plas_carry / count)

#analysis
#focus on T2
phpl_T2 <- filter(phpl, time == "T2")

#random effect of treatment replicate
phpl_T2$rand <- interaction(phpl_T2$evo_line, phpl_T2$rep)

mod1 <- lmer(logden ~ sp * phage * plasmid + (1|rand), data = phpl_T2)
plot(mod1)
mod2 <- lmer(logden ~ sp + phage + plasmid + sp:phage + sp:plasmid + phage:plasmid + (1|rand), data = phpl_T2)
anova(mod1,mod2)
mod3 <- lmer(logden ~ sp + phage + plasmid + sp:plasmid + phage:plasmid + (1|rand), data = phpl_T2)
anova(mod2,mod3) #sig interaction between sp:phage. 0.007785
mod4 <- lmer(logden ~ sp + phage + plasmid + sp:phage + phage:plasmid + (1|rand), data = phpl_T2)
anova(mod2,mod4) #sig interaction between sp:plasmid 0.0201
mod5 <- lmer(logden ~ sp + phage + plasmid + sp:phage + sp:plasmid + (1|rand), data = phpl_T2)
anova(mod2,mod5) #no interaction between phage and plasmid - drop 0.1335
mod6 <- lmer(logden ~ sp + phage + plasmid + sp:plasmid + (1|rand), data = phpl_T2)
anova(mod6,mod5) 
mod7 <- lmer(logden ~ sp + phage + plasmid + sp:phage + (1|rand), data = phpl_T2)
anova(mod7,mod5) 
emmeans::emmeans(mod5, pairwise ~ phage+plasmid|sp) #only sig contrast is between P and PL for O

phpl_T2$evo_line[phpl_T2$evo_line == "PL"] <- "D"

phpl_T2$sp2 <- phpl_T2$sp
phpl_T2$sp2[phpl_T2$sp2 == "P"] <- "(b) P"
phpl_T2$sp2[phpl_T2$sp2 == "O"] <- "(a) O"
phpl_T2$sp2[phpl_T2$sp2 == "V"] <- "(c) V"

bact_denp <- ggplot(phpl_T2, aes(x = evo_line, y = logden)) +
  geom_boxplot(aes(col = evo_line, fill = evo_line)) +
  geom_point(position = position_dodge(0.7), shape = 21, fill = "white", size = 2) +
  facet_wrap(~sp2, ncol = 1) +
  palettetown::scale_fill_poke(pokemon = 'cyndaquil', spread = 6) +
  palettetown::scale_color_poke(pokemon = 'cyndaquil', spread = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "none", strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 12)) +
  scale_x_discrete(labels = c("No Plasmid\nNo Phage", "Plasmid\nNo Phage", "No Plasmid\nPhage", "Plasmid\nPhage")) +
  ylab("Cell density log(CFU/mL)") +
  xlab("Treatment")

library(flextable)
ems <- as.data.frame(emmeans::emmeans(mod5, pairwise ~ phage+plasmid|sp)$contrasts)
ems <- ems[,c(-5)]

ems <- mutate(ems,
                 estimate = round(estimate, digits = 3),
                 SE = round(SE, digits = 3),
                 t.ratio = round(t.ratio, digits = 3),
                 p.value = round(p.value, digits = 3))

emen_tab <- flextable(ems) %>%
  set_header_labels(contrast = "Contrast", sp = "Species", prob = "Proportion", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  hline(i = c(6,12,18), border = officer::fp_border(color="black")) %>%
  align(j = c(1:6), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#relative species proportions
plasdat <- read.csv("github/time_series_week_12.csv", header = T)
plasdat <- filter(plasdat, time == "T2")

#just to grab out totals
cphpl <- read.csv("github/community_structure_090321.csv", header = T)

cphpl_e <- mutate(cphpl,
                  total = P + V + O,
                  P_prop = P / total,
                  V_prop = V / total,
                  O_prop = O / total,
                  lPprop = P_prop * log(P_prop),
                  lVprop = V_prop * log(V_prop),
                  lOprop = O_prop * log(O_prop),
                  lsum = lPprop + lVprop + lOprop,
                  Shann = lsum * -1,
                  even = Shann / log(3))

cphpl_e <- cphpl_e[c(19:24, 13:18, 7:12, 1:6),]
totals <- rep(cphpl_e$total, each = 3)
plasdat$total <- totals

plasdat <- mutate(plasdat,
                  prop = count / total)

plasdat$evo_line[plasdat$evo_line == "P"] <- "PLN"

plasdat$facet <- "(d)"

prop_plot <- ggplot(plasdat, aes(x = evo_line, y = prop, group = interaction(evo_line, sp), col = sp)) +
  geom_boxplot() +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "black", size = 1.5, col = "black") +
  facet_wrap(~facet) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12), strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 12)) +
  palettetown::scale_color_poke(pokemon = "bulbasaur", spread = 3) +
  labs(col = "Species") +
  ylab("Relative proportion") +
  xlab("Treatment") +
  scale_x_discrete(labels = c("No Plasmid\nNo Phage", "Plasmid\nNo Phage", "No Plasmid\nPhage", "Plasmid\nPhage"))+ 
  facet_wrap(~facet)

bact_denp + prop_plot + plot_layout(widths = c(1.5,2))

#analysis of proportions - generalised linear mixed effects model, binomial
y <- cbind(plasdat$count, plasdat$total - plasdat$count)
plasdat$rep2 <- interaction(plasdat$evo_line, plasdat$rep)

mod1 <- glmer(y ~ sp * phage * plasmid + (1|rep2), data = plasdat, family = binomial)
blmeco::dispersion_glmer(mod1) #overdispersion

plasdat <- mutate(plasdat,
                  obs = 1:n())

mod1 <- glmer(y ~ sp * phage * plasmid + (1|rep2) + (1|obs), data = plasdat, family = binomial)
blmeco::dispersion_glmer(mod1)
plot(mod1)
mod2 <- glmer(y ~ sp + phage + plasmid + sp:phage + sp:plasmid + phage:plasmid + (1|rep2) + (1|obs), data = plasdat, family = binomial)
anova(mod1, mod2)
mod3 <- glmer(y ~ sp + phage + plasmid + sp:plasmid + phage:plasmid + (1|rep2) + (1|obs), data = plasdat, family = binomial)
anova(mod2, mod3) #sp:phage interaction
mod4 <- glmer(y ~ sp + phage + plasmid + sp:phage + phage:plasmid + (1|rep2) + (1|obs), data = plasdat, family = binomial)
anova(mod2, mod4) #sp:plasmid interaction
mod5 <- glmer(y ~ sp + phage + plasmid + sp:phage + sp:plasmid + (1|rep2) + (1|obs), data = plasdat, family = binomial)
anova(mod2, mod5) #no phage:plasmid interaction. Drop
mod6 <- glmer(y ~ sp + phage + plasmid + sp:plasmid + (1|rep2) + (1|obs), data = plasdat, family = binomial)
anova(mod5, mod6) #sp:phage
mod7 <- glmer(y ~ sp + phage + plasmid + sp:phage + (1|rep2) + (1|obs), data = plasdat, family = binomial)
anova(mod5, mod7) #sp:plasmid
#model cannot be simplified further
emmeans::emmeans(mod5, pairwise ~ plasmid|sp, type = "response") #for O, proportion is greater with the plasmid. For V, proportion is lower with the plasmid. No effect for P
emmeans::emmeans(mod5, pairwise ~ phage|sp, type = "response") #for O, density lower with phage. For P, density higher with phage
emmeans::emmeans(mod5, pairwise ~ sp|plasmid+phage, type = "response")
emmeans::emmeans(mod5, pairwise ~ sp|plasmid+phage)

#Table for comparing community hierarchy
emens <- emmeans::emmeans(mod5, pairwise ~ sp|plasmid+phage, type = "response")
emens <- as.data.frame(emens)

library(flextable)

#Table S2 - species proportions
emens2 <- emens[c(1:12),-c(4,7)]

emens2 <- mutate(emens2,
                 prob = round(prob, digits = 3),
                 SE = round(SE, digits = 3),
                 asymp.LCL = round(asymp.LCL, digits = 3),
                 asymp.UCL = round(asymp.UCL, digits = 3))

emen_tab <- flextable(emens2) %>%
  set_header_labels(plasmid = "Plasmid", phage = "Phage", sp = "Species", prob = "Proportion", asymp.LCL = "95%LCI", asymp.UCL = "95%UCI") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  hline(i = c(3,6,9), border = officer::fp_border(color="black")) %>%
  align(j = c(1:7), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#Table S3 - pairwise contracts

emens <- emmeans::emmeans(mod5, pairwise ~ sp|plasmid+phage)$contrasts
emens <- data.frame(emens)

emens3 <- emens[,-c(6)]
emens3 <- emens3[,c(2,3,1,4:7)]

emens3 <- mutate(emens3,
                 estimate = round(estimate, digits = 3),
                 SE = round(SE, digits = 3),
                 z.ratio = round(z.ratio, digits = 3),
                 p.value = round(p.value, digits = 3))

emens3$p.value[emens3$p.value == 0] <- '<0.001'

emen_tab2 <- flextable(emens3) %>%
  set_header_labels(phage = "Phage", plasmid = 'Plasmid', contrast = "Contrast", estimate = 'Estimate', z.ratio = "z-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  hline(i = c(3,6,9), border = officer::fp_border(color="black")) %>%
  align(j = c(1:7), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#species evenness in 3-species polycultures

c_even <- read.csv("data/trophic_c_even.csv", header = T)

even1 <- lm(even ~ phage * plasmid, data = c_even)
plot(even1)
even2 <- lm(even ~ phage + plasmid, data = c_even)
anova(even1, even2, test = "F")
even3 <- lm(even ~ plasmid, data = c_even)
anova(even2, even3, test = "F")
even4 <- lm(even ~ phage, data = c_even)
anova(even2, even4, test = "F")
even5 <- lm(even ~ 1, data = c_even)
anova(even4, even5, test = "F")

plotz <- ggplot(data = c_even, aes(x = evo_line, y = even, fill = evo_line, col = evo_line)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  palettetown::scale_fill_poke(pokemon = 'cyndaquil', spread = 6) +
  palettetown::scale_color_poke(pokemon = 'cyndaquil', spread = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "none", strip.background = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  scale_y_continuous(limits = c(0.5, 1.05)) +
  ylab("Species evenness") +
  xlab("Treatment") +
  scale_x_discrete(labels = c("No Plasmid\nNo Phage", "Plasmid\nNo Phage", "No Plasmid\nPhage", "Plasmid\nPhage")) +
  labs(title = "(a)")

#analyse community productivity
prod <- read.csv("data/comms_productivity.csv", header = T)

prod1 <- lm(ltotal ~ phage * plasmid, data = prod)
plot(prod1)
prod2 <- lm(ltotal ~ phage + plasmid, data = prod)
anova(prod1, prod2, test = "F")
prod3 <- lm(ltotal ~ plasmid, data = prod)
anova(prod2, prod3, test = "F")
prod4 <- lm(ltotal ~ phage, data = prod)
anova(prod2, prod4, test = "F")
prod5 <- lm(ltotal ~ 1, data = prod)
anova(prod3, prod5, test = "F")

plotz2 <- ggplot(data = prod, aes(x = evo_line, y = ltotal, fill = evo_line, col = evo_line)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  palettetown::scale_fill_poke(pokemon = 'cyndaquil', spread = 6) +
  palettetown::scale_color_poke(pokemon = 'cyndaquil', spread = 6) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "none", strip.background = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  ylab("Community productivity") +
  xlab("Treatment") +
  scale_x_discrete(labels = c("No Plasmid\nNo Phage", "Plasmid\nNo Phage", "No Plasmid\nPhage", "Plasmid\nPhage")) +
  labs(title = "(b)")

plotz + plotz2 + plot_layout()

#phage and plasmid densities

plas_cs <- read.csv("Plasmid_data.csv", header = T)
plas_cs <- filter(plas_cs, plasmid == "Y")
plas_cs <- mutate(plas_cs,
                  wt = count - plas_carry,
                  prop = plas_carry / count)

plas_cpl <- ggplot(plas_cs, aes(x = sp, y = prop, group = interaction(sp, phage), col = phage, fill = phage))+
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = 'black'), axis.title = element_text(size = 16, colour = 'black'), legend.position = "bottom", legend.title = element_text(size = 14), plot.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_fill_manual(values = c("#274060", "#886020"), labels = c("Plasmid", "Phage & Plasmid")) +
  scale_colour_manual(values = c("#274060", "#886020"), labels = c("Plasmid", "Phage & Plasmid")) +
  ylab("Plasmid proportion") + 
  xlab("Species") +
  labs(col = "Trophic", fill = "Trophic", title = "(b)") 

#model will not fit to standard binomial regression due to high number of 0s and 1s, need to adjust values to help model fitting:
plas_cs$plas_carry2 <- plas_cs$plas_carry
plas_cs$plas_carry2[plas_cs$plas_carry2 == 0] <- 0.001


plas_cs$rep2 <- interaction(plas_cs$evo_line, plas_cs$rep)

mod1 <- glmer(cbind(plas_carry2, count - plas_carry2) ~ sp * phage + (1|rep2), family = binomial, data = plas_cs)
summary(mod1)
blmeco::dispersion_glmer(mod1)
mod2 <- glmer(cbind(plas_carry2, count - plas_carry2) ~ sp + phage + (1|rep2), family = binomial, data = plas_cs)
anova(mod1, mod2)
mod3 <- glmer(cbind(plas_carry2, count - plas_carry2) ~ phage + (1|rep2), family = binomial, data = plas_cs)
anova(mod2, mod3)
mod4 <- glmer(cbind(plas_carry2, count - plas_carry2) ~ sp + (1|rep2), family = binomial, data = plas_cs)
anova(mod2, mod4)
emmeans::emmeans(mod2, pairwise ~ phage, type = "response")
emmeans::emmeans(mod2, pairwise ~ sp, type = "response")

#phage
plas_ph <- read.csv("github/phage_data.csv", header = T)

plas_ph <- mutate(plas_ph,
                  den = dilution * dil_fact * count,
                  logpfus = log10(den))
plas_ph$logpfus[plas_ph$logpfus == -Inf] <- 0

phage_cspl <- ggplot(plas_ph, aes(x = species, y = logpfus, col = plasmid, fill = plasmid, group = interaction(plasmid, species)))  +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black')) +
  theme_bw() +
  scale_fill_manual(values = c("#F8E8A0", "#886020"), labels = c("Phage", "Phage & Plasmid")) +
  scale_colour_manual(values = c("#F8E8A0", "#886020"), labels = c("Phage", "Phage & Plasmid")) +
  theme(axis.text = element_text(size = 14, colour = 'black'), axis.title = element_text(size = 16, colour = 'black'), legend.position = "bottom", legend.title = element_text(size = 14), plot.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  ylab("Phage density log10(PFU/mL)") +
  xlab("Species") +
  labs(col = "Trophic", fill = "Trophic", title = "(a)") +
  scale_x_discrete(labels = c("ORM_57", "CHF7MC", "VAC_51"))

plas_ph$rand <- interaction(plas_ph$line, plas_ph$rep)

cpl1 <- lmer(logpfus ~ species * plasmid + (1|rand), data = plas_ph)
plot(cpl1)
qqnorm(plas_ph$logpfus)
qqline(plas_ph$logpfus, col = "steelblue", lwd = 2)
cpl2 <- lmer(logpfus ~ species + plasmid + (1|rand), data = plas_ph)
anova(cpl1, cpl2)
cpl3 <- lmer(logpfus ~ plasmid + (1|rand), data = plas_ph)
anova(cpl3, cpl2)
cpl4 <- lmer(logpfus ~ species + (1|rand), data = plas_ph)
anova(cpl4, cpl2)
cpl5 <- lmer(logpfus ~ 1 + (1|rand), data = plas_ph)
anova(cpl4, cpl5)

phage_cspl + plas_cpl + plot_layout()


#monoculture data

#analyses looking at effects of phage and plasmid on densities in monoculture

tr_effects <- read.csv("github/monoculture_data.csv", header = T)

tr_plot <- ggplot(tr_effects, aes(x = species, y = int, group = interaction(species, trophic), col = trophic, fill = trophic)) + 
  geom_hline(yintercept=0, linetype = "dashed") +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  scale_fill_manual(values = c("#F8E8A0", "#274060"), labels = c("Phage - No Phage", "Plasmid - No Plasmid")) +
  scale_colour_manual(values = c("#F8E8A0", "#274060"), labels = c("Phage - No Phage", "Plasmid - No Plasmid")) +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = 'black'), axis.title = element_text(size = 16, colour = 'black'), legend.position = "bottom", legend.title = element_text(size = 14), plot.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  ylab("Density change (Trophic - Control)") +
  labs(fill = "Treatment", col = "Treatment") +
  xlab("Species") +
  labs(title = "(a)")

tr_effects$f2 <- interaction(tr_effects$species, tr_effects$trophic)

t_mods <- tr_effects %>%
  nest(-f2) %>%
  mutate(., mod = map(data, ~ t.test(.x$int, mu = 0)))

d_mods_summary <- t_mods %>%
  mutate(params = map(mod, broom::tidy))

d_mods_summary <- unnest(d_mods_summary, params) %>%
  mutate(padjust = p.adjust(p.value, "fdr"))

d_mods_summary$padjust2 <- round(d_mods_summary$padjust, 5) #all species can cope with trophic levels except for P which is significantly negatively affected

library(flextable)
d_mods_summary <- d_mods_summary[,c(1, 5, 6, 8, 9, 13)]

d_mods_summary$Species <- substr(d_mods_summary$f2, 1, 1)
d_mods_summary$f2 <- substr(d_mods_summary$f2, 3, 11)

d_mods_summary <- d_mods_summary[,c(1,7,2,4,5,3,6)]

tab_dat <- d_mods_summary
tab_dat <- mutate(tab_dat,
                  statistic = round(statistic, digits = 3),
                  p.value = round(p.value, digits = 3),
                  padjust2 = round(padjust2, digits = 3),
                  conf.low = round(conf.low, digits = 3),
                  conf.high = round(conf.high, digits = 3))

tab_dat$f2[tab_dat$f2 == "phage"] <- "Phage"
tab_dat$f2[tab_dat$f2 == "plasmid"] <- "Plasmid"

one_flex <- flextable(tab_dat) %>%
  set_header_labels(f2 = "Treatment", statistic = "t-value", p.value = "p-value", padjust2 = "p-value fdr", conf.low = "95%LCI", conf.high = "95%UCI") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() %>%
  hline(i = c(2,4,6), border = officer::fp_border(color="black")) %>%
  align(j = c(1:7), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header')

#phage densities in monoculture

pvo_ph <- read.csv("data/data_interactions/phage_week2_density_PVOphage.csv", header = T)

pvo_ph <- filter(pvo_ph, type == "mono")

phage_plot <- ggplot(pvo_ph, aes(x = species, y = logpfus)) +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = type, col = type)) +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = 'black'), axis.title = element_text(size = 16, colour = 'black'), legend.position = "none", strip.background = element_blank(), plot.title = element_text(size = 16)) +
  ylab("Phage density \nlog10(PFU/mL)") +
  xlab("Species") +
  labs(title = "(b)") +
  scale_fill_manual(values = c("#F8E8A0")) +
  scale_color_manual(values = c("#F8E8A0"))

#plasmid densities in monoculture

pl_m2c <- read.csv("data/plasmid_data.csv", header = T)
pl_m2c <- filter(pl_m2c, ! treatment == "PLP")
pl_m2c[18,5] <- 72

pl_m2c <- filter(pl_m2c, type == "mono")

pl_m2c$prop <- pl_m2c$plas / pl_m2c$total

plasmid_plot <- ggplot(pl_m2c, aes(x = species, y = prop))  +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = type, col = type)) +
  geom_point(position = position_dodge(0.75), shape = 21, fill = "white", size = 2.5, col = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, colour = 'black'), axis.title = element_text(size = 16, colour = 'black'), legend.position = "none", strip.background = element_blank(), plot.title = element_text(size = 16)) +
  labs(col = "Treatment", fill = "Treatment") +
  xlab("Species") + 
  ylab("Plasmid proportion") +
  labs(title = "(c)") +
  scale_fill_manual(values = c("#274060")) +
  scale_color_manual(values = c("#274060"))

pl_m2c$prop[pl_m2c$species == "P"]

mean(c(0.05263158, 0.02777778))

library(patchwork)

antag <- phage_plot + plasmid_plot + plot_layout()

tr_plot + antag + plot_layout(nrow = 2)

#Resistance levels in phage and plasmid-phage communities for Ochrobactrum

res <- read.csv("data/plas_phag_res_comms.csv",header=T)

m1 <- glm(cbind(res,sus) ~ treat, family = binomial, data = res)
summary(m1)
m1 <- glm(cbind(res,sus) ~ treat, family = quasibinomial, data = res)
m2 <- glm(cbind(res,sus) ~ 1, family = quasibinomial, data = res)
anova(m1,m2,test="F")
emmeans::emmeans(m1,pairwise~treat,type="response")

#resistance in mono to polycultures for Ochrobactrum

anc_reso <- read.csv("github/mono_v_comm_ores.csv", header = T)

m1 <- glm(cbind(rs,susc) ~ cond, family = binomial, data = anc_reso)
summary(m1)
m1 <- glm(cbind(rs,susc) ~ cond, family = quasibinomial, data = anc_reso)
m2 <- glm(cbind(rs,susc) ~ 1, family = quasibinomial, data = anc_reso)
anova(m1,m2,test="F")

