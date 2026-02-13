library(ggplot2)
library(data.table)

# Cohen's d = 1
effect_size_target <- 1

# Create grid
d_grid <- data.table(expand.grid(n = seq(1,30,0.01),
                       d = seq(0.001,2,0.001)))

# Compute power for every combination of n and d
d_grid[, power := power.t.test(n = n, d = d, sig.level = 0.05, power = NULL, alternative = "two.sided")$power]
d_1_point <- d_grid[d==effect_size_target&power>0.8][d_grid[d==effect_size_target&power>0.8, which.min(abs(power-0.8)),]]

# Plot power as a function of n and d
ggplot(d_grid, aes(x=n, y=d, fill=power)) +
  geom_raster() +
  scale_fill_viridis_c(begin = 0, end =1, name="Power") +
  scale_x_continuous(breaks=seq(0,30,5), name="Tieranzahl") +
  scale_y_continuous(breaks=seq(0,2,0.2), name="Cohen's d") +
  coord_cartesian(xlim=c(1,30), ylim=c(0,2), expand = F) +
  geom_contour(aes(z=power), color="white") +
  geom_point(data=d_1_point,
             aes(x=n, y=d), color="black", size=1) +
  geom_vline(xintercept = d_1_point$n, linetype="dashed", color="black") +
  geom_hline(yintercept = d_1_point$d, linetype="dashed", color="black") +
  geom_text(aes(x=n+1, y=d+0.1, label=paste0(ceiling(n), " Tiere\nfür Cohen's d = ",effect_size_target,"\nund Power von 80%")),  
            data=d_1_point, hjust="outward", vjust="outward") +
  theme_classic()

ggsave("results/power_B1.png", width=6, height=4)

power.t.test(d = effect_size_target, sig.level = 0.05, power = 0.8, alternative = "two.sided")
