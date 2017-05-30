library(ggplot2)
par(mfrow = c(1,1))

geweke.diag(mcmc(Exponential$gamma))
geweke.diag(mcmc(Exponential$beta))

geweke.diag(mcmc(Weibull$gamma))
geweke.diag(mcmc(Weibull$beta))
geweke.diag(mcmc(Weibull$lambda))

df <- data.frame(x =c('Forest','GDPpc','Deaths'),
                 F =(1/(1+exp(-summary(mcmc(Exponential$gamma))$statistics[-1,1]))),
                 L =(1/(1+exp(-summary(mcmc(Exponential$gamma))$quantiles[-1,1]))),
                 U =(1/(1+exp(-summary(mcmc(Exponential$gamma))$quantiles[-1,5]))))

require(ggplot2)
p <- ggplot(df, aes(x=x, y=F, ymin=L, ymax=U)) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black")) +
  geom_linerange(size=0.5, colour="red") +
  geom_point(size=1, shape=21, fill="black", colour = "white", stroke = 1) +
  geom_hline(aes(yintercept=0), lty=1) +
  coord_flip() +
  ggtitle("Gamma Estimates of z.Exponential Model") +
  labs(y="Change in Prob.", x="")
  
p

df_h <- data.frame(x =c('Forest','GDPpc','Deaths'),
                 F =(exp(summary(mcmc(Exponential$betas))$statistics[-1,1])),
                 L =(exp(summary(mcmc(Exponential$betas))$quantiles[-1,1])),
                 U =(exp(summary(mcmc(Exponential$betas))$quantiles[-1,5])))

h <- ggplot(df_h, aes(x=x, y=F, ymin=L, ymax=U)) +
  theme(axis.text.y   = element_text(size=10),
        axis.text.x   = element_text(size=10),
        axis.title.y  = element_text(size=10),
        axis.title.x  = element_text(size=10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA, color="black")) +
  geom_linerange(size=0.5, colour="red") +
  geom_point(size=1, shape=21, fill="black", colour = "white", stroke = 1) +
  geom_hline(aes(yintercept=1), lty=1) +
  coord_flip() +
  ggtitle("Beta Estimates of z.Exponential Model") +
  labs(y="Hazard Ratio of Civil War Termination", x="")

h
