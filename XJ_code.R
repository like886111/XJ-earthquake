#################################
## IMPORT FUNCTIONS & PACKAGES ##
#################################
library(readr)
library(viridis)
library(R.filesets)
source('/Users/lk/desktop/XJ_functions.R')
source('/Users/lk/desktop/XJ_utilities.R')

#################
## IMPORT DATA ##
#################

# set magnitude of completeness and time domain
#设置完成和时域的大小
M0 = 3.00
my_data <- read.csv("/Users/lk/desktop/xjdz/ht7.csv" , header=TRUE) %>%
  mutate(time_date = as.POSIXct(time_date, format = "%d/%m/%Y %H:%M")) %>% # tranform time column# filer for magnitude
  distinct(time_date, Lat, Lon, .keep_all = TRUE) # remove duplicated lines
# there were 2 duplicated lines
# set starting date (a second before the first event in the sequence)
start.date <- min(my_data$time_date) - 1 

my_data <- my_data %>%
  mutate(time.diff = as.numeric(difftime(time_date, start.date, units = 'days')))

# set up time interval
#设置时间域
T1 = 0
T2 = as.numeric(difftime(max(my_data$time_date) + 1, start.date, units = 'days')) 
# import data
#dd.ama <- read.csv2(file = 'data_M3.0.csv', header = TRUE, sep = ',') %>%
dd.ama<-read.csv(paste("/Users/lk/desktop/xjdz/ht7.csv",sep = ","))%>%
  mutate(time_date = as.POSIXct(paste0(year,'-',month,'-',day,' ',hr,':',min,':',sec)),
         time.diff = as.numeric(difftime(time_date, min(time_date) - 1, units = 'days')),
         Mw = as.numeric(Mw),
         Lon = as.numeric(Lon),
         Lat = as.numeric(Lat),
         Mw.class = cut(Mw, breaks = c(M0, 5, 7))) %>%
  arrange(time_date)
########################
## MCMC MODEL FITTING ##
########################

# 24 mins for 500000
# 16 mins for 300000
#  8 mins for 150000
time.st <- Sys.time()
# generate MCMC posterior samples
MCMC.ama3 <- sampleETASposterior(ts = dd.ama$time.diff, 
                                 magnitudes = dd.ama$Mw, 
                                 M0 = 3.00, 
                                 T=T2, sims = 10000, burnin = 5000, approx = TRUE)
Sys.time() - time.st
#4.264651 mins
# save 
saveRDS(MCMC.ama3, file = '/Users/lk/desktop/MC/MCMC.ht.samples.Rds')
MCMC.ama <- loadRDS('/Users/lk/desktop/MC/MCMC.ht.samples.Rds')

#####
mean(MCMC.ama[,1])
mean(MCMC.ama[,2])
mean(MCMC.ama[,3])
mean(MCMC.ama[,4])
mean(MCMC.ama[,5])

##################################
## MCMC POSTERIOR DISTRIBUTIONS ##
##################################

# parameters name
par.names <- c('mu', 'K', 'alpha', 'c', 'p')


# extrat information for plotting purposes
post.MCMC <- data.frame(value = c(MCMC.ama[,1], MCMC.ama[,2], MCMC.ama[,3], MCMC.ama[,4], MCMC.ama[,5]),
                        param = rep(par.names, each = nrow(MCMC.ama)),
                        type = 'MCMC - posterior')



# data.frame containing priors
prior.MCMC <- rbind(data.frame(x = seq(min(MCMC.ama[,1]), max(MCMC.ama[,1]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dgamma(x, 0.1, 0.1),
                             param = 'mu'),
                    data.frame(x = seq(min(MCMC.ama[,2]), max(MCMC.ama[,2]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'K'),
                    data.frame(x = seq(min(MCMC.ama[,3]), max(MCMC.ama[,3]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'alpha'),
                    data.frame(x = seq(min(MCMC.ama[,4]), max(MCMC.ama[,4]),
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 0, 10),
                             param = 'c'),
                    data.frame(x = seq(min(MCMC.ama[,5]), max(MCMC.ama[,5]), 
                                       length.out = 100)) %>%
                      mutate(pdf = dunif(x, 1, 10),
                             param = 'p')
) %>%
  mutate(type = 'MCMC - prior')

# FIGURE 3
#先验：geom_line(data = prior.MCMC, mapping = aes(x, pdf, color = type,
  #                                         linetype = type) + 
pdf(file = '/Users/lk/desktop/shijianETAS/figure/MCMC.pdf')
ggplot(post.MCMC, aes(x = value)) + geom_density(aes(color = type,
                                                     linetype = type)) + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'bottom')
dev.off()
###########################
## INLABRU MODEL FITTING ##
###########################
# data for Inlabru
data.bru <- data.frame(ts = dd.ama$time.diff, 
                       magnitudes = dd.ama$Mw) %>%
  mutate(idx.p = 1:nrow(dd.ama))

#########copula变换##########
# gamma copula transformation
gamma.t <- function(x, a, b){
  bru_forward_transformation(qgamma, x, a, b)
}

# uniform copula transformation
unif.t <- function(x, a, b){
  bru_forward_transformation(qunif, x, min = a, max = b)
}
# log-gaussian copula transformation
loggaus.t <- function(x, m, s){
  bru_forward_transformation(qlnorm, x, meanlog = m, sdlog = s)
}
# set link functions for Inlabru 
link.f.be <- list(mu = \(x) unif.t(x, 0, 5), 
                  K = \(x) unif.t(x, 0, 5), 
                  alpha = \(x) unif.t(x, 1, 5), 
                  c_ = \(x) unif.t(x, 0,5), 
                  p = \(x) unif.t(x, 1,5))
#
th.init <- list(th.mu = 0,
                th.K = 0,
                th.alpha = 1,
                th.c = 0,
                th.p = 1) 
# options for inlabru 
#inlabru的选项
bru.opt.list <- list(bru_verbose = 4, # type of visual output 
                     bru_max_iter = 100, # maximum number of iterations
                     #bru_method = list(max_step = 0.5),
                     inla.mode = 'experimental', # type of inla algorithm
                     bru_initial = th.init) 

fit_e <- Hawkes.bru(sample.s = data.bru, # data 
                    M0 = M0, # magnitude of completeness
                    T1 = 0, T2 = T2, # time domain
                    link.functions = link.f.be, # link functions
                    coef.t. = 3, # binning parameter (delta)
                    delta.t. = 0.1, # binning parameter (Delta)
                    N.max. = 10, # binning parameter (n.max)
                    bru.opt = bru.opt.list) 
saveRDS(fit_e, file = paste0('/Users/lk/desktop/nihe/ht/fit_e.Rds') )
load("/Users/lk/desktop/nihe/ht/be.Rds")
######################
######拟合优度########
######################
# extract 10000 values from the posterior of the parameters
# this is done in 10 batches of 1000 
# Inlabru replicate case
bru.sample.be <- foreach(i = 1:10, .combine = rbind) %do% {
  sample.p <- generate(fit_e, data.frame(1),  ~ data.frame(mu = link.f.be$mu(th.mu),
                                                            K = link.f.be$K(th.K),
                                                            alpha = link.f.be$alpha(th.alpha),
                                                            c = link.f.be$c_(th.c),
                                                            p = link.f.be$p(th.p)), 
                       n.samples = 1000)
  
  bind_rows(sample.p)
}
saveRDS(bru.sample.be, file = '/Users/lk/desktop/nihe/ht/sample/bru.sample.be.Rds')
bru.be.sample <- loadRDS('/Users/lk/desktop/nihe/ht/sample/bru.sample.be.Rds')

# take last 10000 posterior samples
mcmc.sample.p <- tail(MCMC.ama, 10000)
#初始化空矩阵-元素为Lambda(th)
#不同行对应不同的后验样本
#不同的列对应不同的观察值th
# initialize empty matrices - elements would be Lambda(th)  
# different rows correspond to different posterior samples
# different columns correspond to different observations th
expect.mcmc <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))
expect.bru <- matrix(NA, nrow = 10000, ncol = length(dd.ama$time.diff))
# for each posterior sample
for(i in 1:nrow(bru.sample.be)){
  print(i/nrow(bru.sample.be))
  # find value of Lambda(th) for each time
  expect.mcmc[i,] <- sapply(dd.ama$time.diff, \(x) mcmc.sample.p[i,1]*(x) + sum(exp(
    log.Lambda_h(th = mcmc.sample.p[i,], ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                 mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )
  
  bru.p <- as.numeric(bru.sample.be[i,])
  expect.bru[i,] <- sapply(dd.ama$time.diff, \(x) bru.p[1]*(x) + sum(exp(
    log.Lambda_h2(th = bru.p, ti = dd.ama$time.diff[dd.ama$time.diff <= x],
                  mi = dd.ama$Mw[dd.ama$time.diff <= x], M0 = M0,T1 = 0, T2 = x)))
  )  
}
saveRDS(expect.mcmc, file = '/Users/lk/desktop/nihe/ht/expect/expect.mcmc.Rds')
saveRDS(expect.bru, file = '/Users/lk/desktop/nihe/ht/expect/expect.bru.Rds')

expect.mcmc <- loadRDS('/Users/lk/desktop/nihe/ht/expect/expect.mcmc.Rds')
expect.bru <- loadRDS('/Users/lk/desktop/nihe/ht/expect/expect.bru.Rds')
# 提取Lambda(th)的中位数和分位数值
# extract median and quantiles value of Lambda(th) 
df.res <- data.frame(days = rep(dd.ama$time.diff,2),
                     Lambda.med = c(apply(expect.mcmc, 2, median), 
                                    apply(expect.bru, 2, median)),
                     model = rep(c('MCMC', 'Inlabru'), 
                                 each = nrow(dd.ama)),
                     cumfreq = rep(sapply(dd.ama$time.diff, \(x) sum(dd.ama$time.diff < x)), 2))

# 查找Lambda(th)的累积频率
# find cumulative frequencies for Lambda(th)
xx <- seq(0,1000,length.out = 500)
cs.mcmc <- matrix(NA, ncol = length(xx), nrow = nrow(expect.mcmc))
cs.bru <-  matrix(NA, ncol = length(xx), nrow = nrow(expect.bru))
df.res <- data.frame(days = rep(dd.ama$time.diff,2),
                     Lambda.med = c(apply(expect.mcmc, 2, median), 
                                    apply(expect.bru, 2, median)),
                     Lambda.low = c( apply(expect.mcmc, 2, \(x) quantile(x, 0.025)), 
                                     apply(expect.bru, 2, \(x) quantile(x, 0.025))),
                     Lambda.up = c(apply(expect.mcmc, 2, \(x) quantile(x, 0.975)), 
                                   apply(expect.bru, 2, \(x) quantile(x, 0.975))),
                     model = rep(c('MCMC', 'Inlabru'), 
                                 each = nrow(dd.ama)),
                     cumfreq = rep(sapply(dd.ama$time.diff, \(x) sum(dd.ama$time.diff < x)), 2))
# 查找Lambda(th)的累积频率
# find cumulative frequencies for Lambda(th)
xx <- seq(0,500,length.out = 500)
cs.mcmc <- matrix(NA, ncol = length(xx), nrow = nrow(expect.mcmc))
cs.bru <-  matrix(NA, ncol = length(xx), nrow = nrow(expect.bru))
for(i in 1:nrow(cs.mcmc)){
  cs.mcmc[i,] <- sapply(xx, \(x) sum(expect.mcmc[i,] <= x))
  cs.bru[i,] <- sapply(xx, \(x) sum(expect.bru[i,] <= x))
}
df.cs <- data.frame(Lambda = xx, 
                    cs.med = c( apply(cs.mcmc, 2, median), 
                                apply(cs.bru, 2, median)),
                    cs.low = c( apply(cs.mcmc, 2, \(x) quantile(x, 0.025)), 
                                apply(cs.bru, 2, \(x) quantile(x, 0.025))),
                    cs.up = c( apply(cs.mcmc, 2, \(x) quantile(x, 0.975)), 
                               apply(cs.bru, 2, \(x) quantile(x, 0.975))),
                    model = rep(c('MCMC', 'Inlabru'), 
                                each = ncol(cs.mcmc))
)
# put together
df.total <- data.frame(x = c( df.res$days, df.cs$Lambda ),
                       cumfreq = c( df.res$cumfreq, rep(NA, nrow(df.cs)) ),
                       xx = c( rep(NA, nrow(df.res)), df.cs$Lambda ),
                       med = c( df.res$Lambda.med, df.cs$cs.med ),
                       model = c( df.res$model, df.cs$model ),
                       plot = c( rep('days', nrow(df.res)), rep('Lambda', nrow(df.cs)) )
)


df.total <- data.frame(x = c( df.res$days, df.cs$Lambda ),
                       cumfreq = c( df.res$cumfreq, rep(NA, nrow(df.cs)) ),
                       xx = c( rep(NA, nrow(df.res)), df.cs$Lambda ),
                       med = c( df.res$Lambda.med, df.cs$cs.med ),
                       low = c( df.res$Lambda.low, df.cs$cs.low ),
                       up = c( df.res$Lambda.up, df.cs$cs.up ),
                       model = c( df.res$model, df.cs$model ),
                       plot = c( rep('days', nrow(df.res)), rep('Lambda', nrow(df.cs)) )
)
# FIGURE 2 - TOP ROW
plot.gof.bru <- ggplot(df.total[df.total$model != 'Inlabru',], aes(x = x)) + 
  geom_point(aes(y = cumfreq), size = 1) +
  geom_line(aes(y = xx), linetype = 2) + 
  geom_line(aes(y = med, linetype = model, color = model)) +
  geom_ribbon(aes(ymin = low, ymax = up, 
                  linetype = model, color = model, fill = model), alpha = 0.2) +
  ylab('Cumulative frequencies') + 
  facet_wrap(facets = vars(plot), scales = 'free',
             strip.position = "bottom", 
             labeller = as_labeller(c(days = "days", Lambda = "Lambda") )) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  xlab(NULL) +
  geom_text(data = data.frame(x = c(3000,500), y = 600, plot = c('days', 'Lambda'),
                              lab = c('(a)', '(b)')),
            mapping = aes(x, y, label = lab)) + 
  scale_color_manual(values = gg_color_hue(3)[c(3, 1)]) +
  scale_fill_manual(values = gg_color_hue(3)[c(3, 1)]) + 
  #scale_linetype_manual(values = c(1,2)) + 
  theme(legend.title = element_blank())

# FIGURE 2 - BOTTOM ROW
plot.gof.mcmc <- ggplot(df.total[df.total$model != 'MCMC',], aes(x = x)) + 
  geom_point(aes(y = cumfreq), size = 1) +
  geom_line(aes(y = xx), linetype = 2) + 
  geom_line(aes(y = med, linetype = model, color = model)) +
  geom_ribbon(aes(ymin = low, ymax = up, 
                  linetype = model, color = model, fill = model), alpha = 0.2) +
  ylab('Cumulative frequencies') + 
  facet_wrap(facets = vars(plot), scales = 'free',
             strip.position = "bottom", 
             labeller = as_labeller(c(days = "days", Lambda = "Lambda") )) +
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.placement = "outside") +
  xlab(NULL) +
  geom_text(data = data.frame(x = c(3000,500), y = 600, plot = c('days', 'Lambda'),
                              lab = c('(c)', '(d)')),
            mapping = aes(x, y, label = lab)) + 
  scale_color_manual(values = gg_color_hue(3)[c(2, 1)]) +
  scale_fill_manual(values = gg_color_hue(3)[c(2 ,1)]) + 
  scale_linetype_manual(values = c(1,3)) + 
  theme(legend.title = element_blank())

########################################
#回顾性预测实验#
########################################
# generate posterior sample of ETAS parameters
par.samps <- generate(fit_e, data.frame(), ~ c(th.mu, th.K, th.alpha, th.c, th.p),
                      n.samples = 10000)

# transform them in ETAS scale
par.s.etas <- cbind(mu = link.f.be$mu(par.samps[1,]),
                    K = link.f.be$K(par.samps[2,]),
                    alpha = link.f.be$alpha(par.samps[3,]),
                    c = link.f.be$c_(par.samps[4,]),
                    p = link.f.be$p(par.samps[5,]))
# set number of periods to be forecasts
n.day = 100
# set starting date
f.start <- data.bru$ts[1] + 1e-6
# initialize matrix of forecasting periods
t.lims <- matrix(NA, nrow = n.day, ncol = 3)                                                                                                                                                                                                                                                                                                                                                                                            
# set magnitude above which a forecasting period is splitted
mag.split <- 5.5

# for each period
for(day in 1:n.day){
  # set starting and end time of the forecasting period 
  if(day == 1){
    f.T1 <- f.start
    f.T2 <- f.start + 1
  }
  else{
    f.T1 = f.T2
    f.T2 = f.T1 + 1
  }
  # select data in the forecasting period
  data.T1.T2 <- data.bru[data.bru$ts > f.T1 & data.bru$ts < f.T2, ]
  # check if any event above mag.split
  if(any(data.T1.T2$magnitudes > mag.split)){
    f.T2 = data.T1.T2$ts[which.max(data.T1.T2$magnitudes > mag.split)] + 1e-6
  }
  # store time forecasting interval
  t.lims[day, ] <- c(f.T1,  ( (f.T2 + f.T1)/2 ), f.T2)
  # produce forecast for the period
  single.fore <- cat_forecast(theta.samp = par.s.etas,
                              fore.T1 = f.T1,
                              fore.delta = (f.T2 - f.T1),
                              M0 = M0, beta.p = beta.ml,
                              data.input = data.bru,
                              folder.path = '/Users/lk/desktop/nihe/ht/fore_cat100/',
                              fore.name = paste0('fore.day.', day))
}

# retrieve daily forecasts and get quantiles of the number of events per day
# retrieve daily forecasts and get quantiles of the number of events per day
# retrieve daily forecasts and get quantiles of the number of events per day
N.quant.sim <- foreach(day = 1:n.day, .combine = rbind) %do% {
  print(day)
  f.cat <- read.table(file = paste0('/Users/lk/desktop/nihe/ht/fore_cat100/fore.day.', day ,'.txt'),
                      header = TRUE)
  n.sim <- vapply(1:nrow(par.s.etas), \(x) sum(f.cat$cat.idx == x), 0)
  n.true <- sum(data.bru$ts > t.lims[day, 1] & 
                  data.bru$ts <= t.lims[day,3] )
  data.frame(lower = quantile(n.sim, 0.025),
             median = quantile(n.sim, 0.5),
             upper = quantile(n.sim, 0.975),
             true = n.true)
}

# set rownames
rownames(N.quant.sim) <- NULL
# create data.frame of mids points of each time period for plotting
t.lims.df <- data.frame(t.mid = t.lims[,2])

periods <- 1:n.day
# FIGURE 4 TOP ROW
plot.foreN.natural <- ggplot(N.quant.sim, aes(x = t.lims.df$t.mid, 
                                              y = median)) + 
  geom_line(color = 'red') + 
  geom_ribbon(aes(xmin = t.lims.df$t.mid, 
                  xmax = t.lims.df$t.mid, 
                  ymin = lower, 
                  ymax = upper),
              alpha = 0.2, color = 'orange', fill = 'orange') + 
  geom_point(aes(x = t.lims.df$t.mid, 
                 y = true)) + #, size = 0.5)) + 
  #scale_y_log10() + 
  xlab('periods') + 
  ylab('N') + 
  annotate('text', x = 30, y = 10, label = '(a)')+
  theme_classic()

# FIGURE 4 BOTTOM ROW
# Remove periods with 0 events
idx.good <- N.quant.sim$true > 0 
plot.foreN.log <- 
  ggplot(N.quant.sim[idx.good,], aes(x = t.lims.df$t.mid[idx.good], 
                                     y = log(median))) + 
  geom_line(color = 'red') + 
  geom_ribbon(aes(xmin = t.lims.df$t.mid[idx.good], 
                  xmax = t.lims.df$t.mid[idx.good], 
                  ymin = log(lower), 
                  ymax = log(upper)),
              alpha = 0.2, color = 'orange', fill = 'orange') + 
  geom_point(aes(x = t.lims.df$t.mid[idx.good], 
                 y = log(true))) + #, size = 0.5)) + 
  #scale_y_log10() + 
  xlab('periods') + 
  ylab('log(N)') + 
  annotate('text', x = 85, y = log(10), label = '(b)')+
  theme_classic()
# FIGURE 4
pdf('/Users/lk/desktop/figure00.pdf')
multiplot(plot.foreN.natural, plot.foreN.log)
dev.off()