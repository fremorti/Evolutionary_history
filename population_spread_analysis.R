###################################
#                                 #
#  Analysis of population spread  #
#                                 #
###################################
library(car)
library(tidyverse)
library(ggplot2)
library(brms)
library(gridExtra)

#reading-in spread data
data <- read.csv('data_mesocosms.csv', header = 1)
data$metapop <- as.factor(paste0(data$treatment, '.', data$metarep))
data$treat8 <- data$treatment == 8
data$treat16 <- data$treatment == 16
data$treatment <- as.factor(data$treatment)


#exploratory plots
#spread per treatment
ggplot(data1, aes(x = day, y = edge, color = treatment))+
  geom_point()+
  geom_smooth(method = 'lm')
#spread per mesocosm
ggplot(data1, aes(x = day, y = edge, color = treatment, group = metapop))+
  geom_point()+
  geom_smooth(method = 'lm')

#model spread as a function of the evolutionary history (metapop)
ffedge <- brmsformula(edge  ~ 0 + Intercept + treatment*day + (1 + day|metapop))
fedge <- brm(data = data, family = 'gaussian',
             formula = ffedge,
             prior = c(prior(normal(0,4), class = b),
                       prior(cauchy(0,2), class = sd),
                       prior(lkj(2), class = cor)),
             iter = 5000, warmup = 2000, chains = 2, cores = 2)
fedge$fit

#calculate sd of parameters
mcmc <- as.matrix(fedge)

#sd of day
sd.d <- sd(data1$day)
sd.day <- abs(mcmc[,4])*sd.d

#sd of treatment
data.b <- cbind(data1$treatment == '8', data1$treatment == '16')
sd.b <- apply(mcmc[, 2:3]%*%t(data.b), 1, sd)
bs <- cbind(mcmc[,2:3], rep(0, 6000))

#sd of treatment:day
data.bdayp <- t(data.b*data1$day)
sd.bday <- apply(mcmc[, 5:6] %*% data.bdayp, 1, sd)

#sd of metapop
data.meta <- sapply(data1$metapop, function(x) levels(data1$metapop) == x)
sd.metapop <- apply(mcmc[, 11:25] %*% data.meta, 1, sd)

#sd of metapop:day
data.metadayp <- data.meta*data1$day
sd.metaday <- apply(mcmc[, 26:40] %*% data.metadayp, 1, sd)

#sd of residuals
sd.resid <- apply(residuals(fedge, summary = FALSE), 1, sd)

sd.all <- cbind(sd.b, sd.day, sd.bday, sd.metapop, sd.metaday, sd.resid)%>%as_tibble()%>%rename('con' = sd.b, 'day' = 'sd.day', 'con:day' = 'sd.bday', 'mesoc' = 'sd.metapop', 'mesoc:day' = 'sd.metaday', 'resid' = 'sd.resid')
summary(sd.all)

#plot all sources of variation in the model
sd.all_ <- sd.all%>%pivot_longer(cols = colnames(.))
ggplot(sd.all_, aes(x = name, y = value))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), scale = 'width', color = 'grey40', aes(fill = name))+
  scale_fill_brewer( palette = "Set2")+
  xlab('level of the model')+
  ylab('standard deviation of parameters')+
  ggtitle('sources of variation in population edge')+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))
  


#calculate variaiton proportional to the residual and calculate proportional differences
sd.all$b_ <- sd.b/sd.resid
sd.all$day_ <- sd.day/sd.resid
sd.all$bday_ <- sd.bday/sd.resid
sd.all$metapop_ <- sd.metapop/sd.resid
sd.all$metaday_ <- sd.metaday/sd.resid
sd.all$'b:day_metapop:day' <- sd.bday/sd.metaday
sd.all$'b_metapop' <- sd.metapop/sd.b

#plot proportions involving slopes of treatment and mesocosm
ggplot(sd.all[c(9,11,12)]%>%select('b:day_metapop:day' = 'b:day_metapop:day', 'b:day_resid' = bday_, 'metapop:day_resid' = metaday_)%>%pivot_longer(cols = colnames(.)), aes(name, value))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), scale = 'width', fill = 'bisque3', color = 'bisque4')+
  geom_hline(yintercept = 1, linetype = 'dashed', size = 1, color = 'bisque4')+
  scale_x_discrete(labels = c('sd_conday / sd_mesocday', 'sd_conday / sd_resid', 'sd_mesocday / sd_resid'))+
  xlab('contrasted levels')+
  ylab('proportional difference standard deviation of parameters')+
  ggtitle('proportional differences in sources of variance: slopes')+
  theme_minimal()+
  theme(legend.position = "none")
#plot proportions involving intercepts of treatment and mesocosm
ggplot(sd.all[c(7,10,13)]%>%pivot_longer(c(1, 2, 3))%>%mutate(name = factor(.$name, levels = c('b_', 'metapop_', 'b_metapop'))), aes(name, value))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), scale = 'width')+
  geom_hline(yintercept = 1, linetype = 'dashed', size = 1, color = 'darkred')+
  scale_x_discrete(labels = c('sd_con / sd_resid', 'sd_mesoc / sd_resid', 'sd_mesoc / sd_con'))+
  ylim(c(-0.1, 3))+
  xlab('contrasted levels')+
  ylab('proportion in standard deviation of estimated parameters')+
  ggtitle('proportional differences in sources of variance: intercepts')
  
         



#visualizing the estimated treatment effect (although small-ish estimated variance component)
#culculate group estimates from posterior
pfedge <- fedge %>% posterior_samples() %>% mutate(int_8 = .$b_Intercept + .$b_treatment8,
                                                     int_16 = .$b_Intercept + .$b_treatment16,
                                                     slo_8 = .$b_day + .$'b_treatment8:day',
                                                     slo_16 = .$b_day + .$'b_treatment16:day',
                                                     int_8_16 = .$b_treatment8 - .$b_treatment16,
                                                     slo_8_16 = .$'b_treatment8:day'-.$'b_treatment16:day')

pslo1 <- pfedge %>% select('4' = b_day, '8' = slo_8, '16' = slo_16) %>% pivot_longer(cols = c('4', '8', '16'), 'treatment')%>%
  mutate(treatment = factor(.$treatment, levels = c('4', '8', '16')))
#plot estimated slopes
(fedge1_splot <- ggplot(pslo1, aes(x = treatment, y = value, color = treatment))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
    scale_color_brewer( palette = "Dark2")+
    theme(legend.position='none')+
    ylab('estimated slope')+
    xlab('connectedness treatment'))
pslo1_d <- pfedge %>% select('8_4' = 'b_treatment8:day', '16_4' = 'b_treatment16:day', '8_16' = slo_8_16) %>% pivot_longer(cols = c('8_4', '16_4', '8_16'),'difference')%>%
  mutate(difference = factor(.$difference, levels = c('8_4', '16_4', '8_16')))
#plot estimated differences in slopes
(fedge1_splot_d <- ggplot(pslo1_d, aes(x = difference, y = value))+
    geom_violin(draw_quantiles = c(0.09, 0.5, 0.91))+
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkred', size = 1)+
    theme(legend.position='none')+
    ylab('estimated difference in slope'))

nde <- tibble(day = rep(1:27, 3), treatment = rep(c(4,8,16), each = 27))
postf1 <- fitted(fedge, newdata = nde, re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nde)%>%mutate(treatment = as.factor(.$treatment))

#plot highest likelyhood slopes per treatment (and 91% spread)
(e3 <- ggplot(postf1, aes(x = day, y = Estimate, group = treatment))+
    geom_jitter(data = data, aes(x = day, y = edge, color = treatment), width = 0.4, height = 0.4)+
    geom_line(aes(color = treatment), size = 1.25)+
    geom_ribbon(aes(ymin = Q9, ymax = Q91, fill = treatment), alpha = 0.25) +
    scale_fill_brewer( palette = "Dark2")+
    scale_color_brewer( palette = "Dark2")+
    ylab('population edge')+
    labs(color = 'connectedness')+
    guides(fill = FALSE)+
    ggtitle('connectedness effect on spread')+
    theme(legend.position='top'))


grid.arrange(grobs = list(e3, fedge1_splot, fedge1_splot_d), widths = c(1,1,2), layout_matrix = rbind(c(1,1, 2),c(1,1, 3)))


##################
# trait analysis #
##################

#load trait data
LH <- read.table('LHtraits.txt', header = TRUE)
LH$repr <- LH$JUVEN + LH$AFEMA + LH$AMALES
LH$metapop <- paste0(LH$TREAT, '.', substring(LH$PATCH, nchar(as.character(LH$PATCH)), nchar(as.character(LH$PATCH))))
g <- as_tibble(LH)%>%filter(TIME == 10)%>%select(metapop, repr)%>%group_by(metapop)%>%summarize(repr = mean(repr))
data2 <- merge(data1, g, all.x = TRUE)
disp <- read.table('dispersal.txt', header = TRUE)
disp$metapop <- paste0(disp$TREAT, '.', substring(disp$REPLICA, nchar(as.character(disp$REPLICA)), nchar(as.character(disp$REPLICA))))
d <- as_tibble(disp)%>%group_by(metapop)%>%summarise(end = sum(End), start = mean(Density))%>%mutate(disp = end/start)
data3 <- merge(data2, d[c(1,4)], all.x = TRUE)

#model spread as dependent on reproductive success
ffLH <- brmsformula(edge  ~ 0 + Intercept + repr*day + (1+day|metapop))
get_prior(ffLH, data = data3)
fLH <- brm(data = data3, family = 'gaussian',
           formula = ffLH,
           prior = c(prior(normal(0,4), class = b),
                     prior(cauchy(0,2), class = sd),
                     prior(lkj(2), class = cor)),
           iter = 5000, warmup = 2000, chains = 2, cores = 2)
fLH$fit

#calculate variation components
LHmcmc <- as.matrix(fLH)
LHdata <- droplevels(data3[!is.na(data3$repr),])

LHsd.d <- sd(LHdata$day)
LHsd.day <- abs(LHmcmc[,3])*sd.d

LHsd.r <- sd(LHdata$repr)
LHsd.repr <- abs(LHmcmc[,2])*LHsd.r

LHsd.rd <- sd(LHdata$repr*LHdata$day)
LHsd.reprday <- abs(LHmcmc[,4])*LHsd.rd

LHdata.meta <- sapply(LHdata$metapop, function(x) levels(LHdata$metapop) == x)
LHsd.metapop <- apply(LHmcmc[, 9:21] %*% LHdata.meta, 1, sd)

LHdata.metadayp <- LHdata.meta*LHdata$day
LHsd.metaday <- apply(mcmc[, 22:34] %*% LHdata.metadayp, 1, sd)

LHsd.resid <- apply(residuals(fLH, summary = FALSE), 1, sd)

LHsd.all <- cbind(LHsd.day, LHsd.repr, LHsd.reprday, sd.metapop, sd.metaday, sd.resid)%>%as_tibble()%>%rename('day' = 'LHsd.day', 'repr' = LHsd.repr, 'repr:day' = LHsd.reprday, 'metapop' = 'sd.metapop', 'metapop:day' = 'sd.metaday', 'resid' = 'sd.resid')
summary(LHsd.all)

#plot sources of variation
LHsd.all_ <- LHsd.all%>%pivot_longer(cols = colnames(.))
LHsd.all_$name <- factor(LHsd.all_$name, levels = c('day', 'repr', 'repr:day', 'metapop', 'metapop:day', 'resid'))
ggplot(LHsd.all_, aes(x = name, y = value))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), scale = 'width', aes(fill = name))+
  scale_x_discrete(labels = c('day', 'repr', 'repr:day', 'mesoc', 'mesoc:day', 'resid'))+
  scale_fill_brewer( palette = "Set2")+
  xlab('level of the model')+
  ylab('standard deviation of parameters')+
  ggtitle('sources of variance in population edge')+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))




#calculate estimated slopes and intercepts from posterior
nLH <- tibble(day = rep(1:27, 71), repr = rep(seq(20, 55, length.out = 71), each = 27))
postLH <- fitted(fLH, newdata = nLH, nlpar = "b", re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nLH)%>%filter(day == 1)
#plot estimated intercepts
(LH1 <- ggplot(postLH, aes(x = repr, y = Estimate))+
  geom_line(size = 1.25)+
  geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
  ylab('estimated intercept')+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))+
  ggtitle('reproductive success effect on initial edge'))

postLH_ <- fitted(fLH, newdata = nLH, nlpar = "a", re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nLH)%>%filter(day == 1)
#plot estimated slopes
(LH1_ <- ggplot(postLH_, aes(x = repr, y = Estimate))+
  geom_line(size = 1.25)+
  geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
  ylab('estimated slope')+
  ggtitle('reproductive success effect on spread rate'))+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))


#model population spread as dependent on dispersal capacity
ffd <- brmsformula(edge  ~ 0 + Intercept + disp*day + (1+day|metapop))
get_prior(ffd, data = data3)
fd <- brm(data = data3, family = 'gaussian',
          formula = ffd,
          prior = c(prior(normal(0,4), class = b),
                    prior(cauchy(0,2), class = sd),
                    prior(lkj(2), class = cor)),
          iter = 5000, warmup = 2000, chains = 2, cores = 2)
fd$fit

#calculate components of variation of the model
ffd <- brmsformula(edge  ~ day*disp)
bayes_R2(fd)

dmcmc <- as.matrix(fd)
ddata <- droplevels(data3[!is.na(data3$disp),])

dsd.d <- sd(ddata$day)
dsd.day <- abs(dmcmc[,3])*sd.d

dsd.di <- sd(ddata$disp)
dsd.disp <- abs(dmcmc[,2])*dsd.di

dsd.dd <- sd(ddata$disp*ddata$day)
dsd.dispday <- abs(dmcmc[,4])*dsd.dd

ddata.meta <- sapply(ddata$metapop, function(x) levels(ddata$metapop) == x)
dsd.metapop <- apply(dmcmc[, 9:21] %*% ddata.meta, 1, sd)

ddata.metadayp <- ddata.meta*ddata$day
dsd.metaday <- apply(mcmc[, 22:34] %*% ddata.metadayp, 1, sd)

dsd.resid <- apply(residuals(fd, summary = FALSE), 1, sd)

dsd.all <- cbind(dsd.day, dsd.disp, dsd.dispday, sd.metapop, sd.metaday, sd.resid)%>%as_tibble()%>%rename('day' = 'dsd.day', 'disp' = dsd.disp, 'disp:day' = dsd.dispday, 'metapop' = 'sd.metapop', 'metapop:day' = 'sd.metaday', 'resid' = 'sd.resid')
summary(dsd.all)

#plot sources of variation in the model
dsd.all_ <- dsd.all%>%pivot_longer(cols = colnames(.))
ggplot(dsd.all_, aes(x = name, y = value))+
  geom_violin(draw_quantiles = c(0.09, 0.5, 0.91), scale = 'width', aes(fill = name))+
  scale_x_discrete(labels = c('day', 'disp', 'disp:day', 'mesoc', 'mesoc:day', 'resid'))+
  scale_fill_brewer( palette = "Set2")+
  xlab('level of the model')+
  ylab('standard deviation of parameters')+
  ggtitle('sources of variance in population edge')+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))


nd <- tibble(day = rep(1:27, 50), disp = rep(seq(0.1, 0.6, length.out = 50), each = 27))
postd <- fitted(fd, newdata = nd, nlpar = "b", re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nd)%>%filter(day == 1)
postd_ <- fitted(fd, newdata = nd, nlpar = "a", re_formula = NA, probs = c(0.09, 0.91))%>%
  as_tibble()%>% bind_cols(nd)%>%filter(day == 1)

#plot estimated intercepts
(d1 <- ggplot(postd, aes(x = disp, y = Estimate))+
  geom_line(size = 1.25)+
  geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
  ylab('estimated intercept')+
  ggtitle('dispersal effect on initial edge'))+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))

#plot estimated slopes
(d1_ <- ggplot(postd_, aes(x = disp, y = Estimate))+
  geom_line(size = 1.25)+
  geom_ribbon(aes(ymin = Q9, ymax = Q91), fill = "black", alpha = 0.1) +
  ylab('estimated slope')+
  ggtitle('dispersal effect on spread rate'))+
  theme_minimal()+
  theme(legend.position = "none", text = element_text(size = 15))
