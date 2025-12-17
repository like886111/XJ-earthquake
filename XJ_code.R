library(readr)
library(viridis)
library(R.filesets)
source('/Users/lk/desktop/XJ_functions.R')
source('/Users/lk/desktop/XJ_utilities.R')

M0 = 3.00
my_data <- read.csv('' , header=TRUE) %>%
  mutate(time_date = as.POSIXct(time_date, format = "%d/%m/%Y %H:%M")) %>% # tranform time column# filer for magnitude
  distinct(time_date, Lat, Lon, .keep_all = TRUE) # remove duplicated lines
# there were 2 duplicated lines
# set starting date (a second before the first event in the sequence)
start.date <- min(my_data$time_date) - 1 

my_data <- my_data %>%
  mutate(time.diff = as.numeric(difftime(time_date, start.date, units = 'days')))

T1 = 0
T2 = as.numeric(difftime(max(my_data$time_date) + 1, start.date, units = 'days')) 
# import data
#dd.ama <- read.csv2(file = 'data_M3.0.csv', header = TRUE, sep = ',') %>%
dd.ama<-read.csv(paste('',sep = ","))%>%
  mutate(time_date = as.POSIXct(paste0(year,'-',month,'-',day,' ',hr,':',min,':',sec)),
         time.diff = as.numeric(difftime(time_date, min(time_date) - 1, units = 'days')),
         Mw = as.numeric(Mw),
         Lon = as.numeric(Lon),
         Lat = as.numeric(Lat),
         Mw.class = cut(Mw, breaks = c(M0, 5, 7))) %>%
  arrange(time_date)

data.bru <- data.frame(ts = dd.ama$time.diff, 
                       magnitudes = dd.ama$Mw) %>%
  mutate(idx.p = 1:nrow(dd.ama))

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
v.par <- c(28,29,30,31)

# set link function list 
link.f.gamma.list <- lapply(1:length(v.par), \(i)  
                            list(mu = \(x) gamma.t(x, force(v.par[i]), 
                                                   10*force(v.par[i])), 
                                 K = \(x) gamma.t(x, 2*force(v.par[i]), 
                                                  force(v.par[i])), 
                                 alpha = \(x) gamma.t(x, 2*force(v.par[i]), 
                                                      force(v.par[i])), 
                                 c_ = \(x) gamma.t(x, force(v.par[i]), 
                                                   10*force(v.par[i])), 
                                 p = \(x) 1 + gamma.t(x, force(v.par[i]), 
                                                      2*force(v.par[i])))
)


# set Inlabru options list
bru.opt.list.gamma <- list(bru_verbose = 3,
                           bru_max_iter = 100,
                           #bru_method = list(max_step = 0.5),
                           inla.mode = 'experimental',
                           bru_initial = list(th.mu = 0,
                                              th.K = 0,
                                              th.alpha = 1,
                                              th.c = 0,
                                              th.p = 1))
# fit the model for each parameter value
for(i in 1:length(v.par)){
  # extract link functions
  link_f <- link.f.gamma.list[[i]]
  # extract gamma value
  gamma.par <- v.par[i]
  prior.name <- paste0('gamma par : ', gamma.par)
  cat(prior.name, '\n')
  # fit the model 
  fit_ <- Hawkes.bru(sample.s = data.bru, 
                     M0 = M0,
                     T1 = 0, T2 = T2, 
                     coef.t. = 2, 
                     delta.t. = 0.1, 
                     N.max. = 3, 
                     link.functions = link_f,
                     bru.opt = bru.opt.list.gamma)
  # save
  saveRDS(fit_, file = paste0('', 'fit_par', gamma.par, '_N.max10.Rds') )
  cat('Perc compl', i/length(v.par), '\n')
}

#prior/gamma
# takes file names
l.files <- list.files(path = '')
# load model and take parameter values from file name
fit.prior.list <- foreach(i = 1:length(l.files)) %do% {
  fit_ <- loadRDS(paste0('', l.files[i]))
  g.p <- parse_number(l.files[i])
  list(fit = fit_,
       gamma.par = g.p)
}

# create list of posterior distributions
post.pr.list <- lapply(1:length(fit.prior.list), \(idx) 
                       extract.post.df(fit.prior.list[[idx]]$fit, 
                                       link.f.gamma.list[[idx]]) %>%
                         mutate(gamma.p = fit.prior.list[[idx]]$gamma.par))

# bind rows and set priors
post.pr.bind <- bind_rows(post.pr.list) 
post.pr.bind <- post.pr.bind %>% 
  mutate(prior = case_when(param == 'mu' ~ dgamma(x, gamma.p, 10*gamma.p),
                           param == 'K' ~ dgamma(x, 2*gamma.p, gamma.p),
                           param == 'alpha' ~ dgamma(x, 2*gamma.p, gamma.p),
                           param == 'c' ~ dgamma(x, gamma.p, 10*gamma.p),
                           param == 'p' ~ dgamma(x - 1, gamma.p, 5*gamma.p)))

post.pr.bind <- post.pr.bind %>%
  filter(param %in% c("alpha", "p"))
# FIGURE 7
ggplot(post.pr.bind, 
       aes(x,y, color = factor(gamma.p))) + 
  geom_line(aes(linetype = 'Posterior')) + 
  geom_line(aes(y = prior, linetype = 'Prior')) + 
  labs(color = ~ gamma, linetype = '') + 
  scale_y_log10() +
  ylab('log10(pdf)') + 
  xlab('value') + 
  facet_wrap(facets = vars(param), scales = 'free',
             labeller = label_parsed) +
  theme_classic() + 
  theme(legend.position = 'bottom') + 
  scale_color_viridis(discrete = TRUE)













