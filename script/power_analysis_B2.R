library(lme4)
library(simr)
library(data.table)
library(ggplot2)

# design parameters
p_vals <- function(n) {
  print(n)
  IRR <- 1.63
  KO_lambda <- 120
  sessions <- 5
  design <- expand.grid(
    subject   = rep(factor(1:(2 * n)), times = sessions),
    condition = factor(c("pre", "during"), levels = c("pre", "during"))
  )
  
  # make design table
  set.seed(1) # set seed for reproducibility
  design <- data.table(design,
                       group=factor(c("WT", "KO"), levels=c("WT", "KO")))
  design[,lambda:=0,][
    group=="KO", lambda:=lambda+log(KO_lambda),][
      condition=="during"&group=="KO", lambda:=lambda+log(IRR),]
  design[,count:=rpois(.N, exp(lambda))][,lambda:=NULL]

  # define model
  model <- glmer(
    count ~ group * condition + (1 | subject),
    family = poisson,
    data   = design
  )
  
  # set fixed effects and random effects variance
  fixef(model) <- c(
    "(Intercept)" = 0,                    # WT log(baseline rate)
    "groupKO" = log(KO_lambda),                 # KO log(baseline rate)
    "conditionduring" = 0,                # WT log(IRR)
    "groupKO:conditionduring" = log(IRR) # KO log(IRR)
  )
  
  # variance of random effect
  VarCorr(model)$subject[1] <- 0.4^2
  
  # compare to alternative model and get p-values
  powerOut <- powerSim(
    model,
    test = fcompare(count ~ group + condition, "lr"),
    nsim = 1000,
  )
  return(data.table(n_animals = n, power = powerOut$pval))
}

n_animals <- c(2,5,10,12:16,20)
power_summary <- lapply(n_animals, p_vals)
power_dt <- rbindlist(power_summary)[,.(n = sum(power < 0.05)), by=.(n_animals)]
power_dt[,`:=`(Power = n/1000, CI = list(qbeta(c(0.025, 0.5,0.975), shape1 = 1+n, shape2 = 1000-n))), by=.(n_animals)]
power_dt[,`:=`(CI_low = unlist(CI)[1], CI_mid = unlist(CI)[2], CI_high = unlist(CI)[3]), by=.(n_animals)][,CI:=NULL,]

# save results as csv
fwrite(power_dt, file="results/power_B2.csv")

# plot power curve
ggplot(power_dt, aes(x=n_animals, y=Power)) +
  geom_line()+
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high), alpha=0.2)+
  scale_x_continuous(breaks=seq(0, 20, 5), name = "Tierzahl (N)")+
  scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1), name="Power")+
  geom_hline(yintercept = 0.8, linetype="dashed")+
  geom_vline(xintercept = power_dt[CI_low>0.8,min(n_animals),], linetype="dashed")+
  annotate("text", x=power_dt[CI_low>0.8,min(n_animals),]-1, y=0.1, label=paste("  ",power_dt[CI_low>0.8,min(n_animals),],  "Tiere\nfür 80% Power"), hjust = "inward", vjust="outward")+
  coord_cartesian(ylim = c(0,1), xlim = c(0,20), expand = F)+
  theme_classic()

# save plot
ggsave("results/power_B2.png", width=6, height=4)
