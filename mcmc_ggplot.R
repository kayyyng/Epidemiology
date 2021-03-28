subCladeSize = c("2", "3", "5")

library(ggplot2)
library(gridExtra)

for (i in 1:length(subCladeSize)) {
  df = read.csv(paste("mcmc_result/strict/equal_more_", subCladeSize[i], 
                      "/all_country_mcmc_res.csv", sep = ""), header = FALSE)
  df_2 = read.csv(paste("../Backups/mcmc_result_before_deconv/strict/equal_more_",
                    subCladeSize[i], "/all_country_mcmc_res.csv", sep = ""), header = FALSE)
  colors = c("Without deconvolution" = "black", "With deconvolution" = "red")
  plot = ggplot(df, aes(x = as.numeric(row.names(df)), y = V1)) + 
    geom_line(data = df_2, aes(x = as.numeric(row.names(df)), y = V1, color = "Without deconvolution")) +
    geom_line(aes(group=1, color = "With deconvolution")) + 
    labs(x = "", y = "", title = paste("Subclade Size = ", subCladeSize[i], sep = ""), color = "Legend") + 
    scale_x_continuous(limits = range(0,30000)) + scale_y_continuous(breaks=seq(1.1,1.362,0.05)) +
    scale_color_manual(values = colors)
  name = paste("p", i, sep = "")
  assign(name, plot)
}

grid.arrange(p1, p2, p3, ncol= 1)
