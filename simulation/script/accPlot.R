readResult = function(fold, method, dim, col = 1){
  
  file2 = paste0(fold, "/P_", dim, "/prop_0.2/summary/", method, ".txt")
  x2 = read.table(file2, sep =" ", header = T)
  x2 = c(mean(x2[, col], na.rm = T), sd(x2[, col], na.rm = T))
  
  file3 = paste0(fold, "/P_", dim, "/prop_0.4/summary/", method, ".txt")
  x3 = read.table(file3, sep =" ", header = T)
  x3 = c(mean(x3[, col], na.rm = T), sd(x3[, col], na.rm = T))
  
  file4 = paste0(fold, "/P_", dim, "/prop_0.6/summary/", method, ".txt")
  x4 = read.table(file4, sep =" ", header = T)
  x4 = c(mean(x4[, col], na.rm = T), sd(x4[, col], na.rm = T))
  
  file5 = paste0(fold, "/P_", dim, "/prop_0.8/summary/", method, ".txt")
  x5 = read.table(file5, sep =" ", header = T)
  x5 = c(mean(x5[, col], na.rm = T), sd(x5[, col], na.rm = T))
  
  result = rbind(x2, x3, x4, x5)
  colnames(result) = c("Measure", "SD")
  return(result)
}



plot_acc <-function(folder = "simulation2/patternMOB1", dim = 500, legend_pos = "none",
                    title = TeX("$p = 500, \\Delta = 1$")){
  
  metric1 = readResult(folder, "accuracy", dim, col=1)
  metric2 = readResult(folder, "accuracy", dim, col=2)
  metric3 =  readResult(folder, "accuracy", dim, col=3)

  df = rbind(metric1, metric2, metric3)
  method = factor(rep(c("Sensitivity", "Specificity", "Accuracy"), each = 4), 
                  levels = c("Sensitivity", "Specificity", "Accuracy"))
  df = data.frame(df, ratio=rep((0.2*1:4), 3), method = method)
  
  
  library(ggplot2)
  library(latex2exp)
  
  plotARI <- ggplot(data=df, mapping = aes(x=ratio, y = Measure)) + 
    geom_line(aes(color = method)) + 
    geom_point(aes(color = method)) +
    ylim(-0.05, 1.05) + 
    geom_errorbar(aes(ymin = Measure-SD, ymax = Measure + SD, color = method), width = 0.01)+
    ggtitle(title) + theme_bw() + 
    theme(legend.position = legend_pos, 
          legend.spacing = unit(0.1, "cm"),
          legend.key.height = unit(10, "pt"),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
          axis.title.x = element_blank(), 
          axis.title.y = element_text(face = "bold")) 
  
  plotARI
  
}

p1 = plot_acc(folder = "simulation3/patternMOB1", dim = 500, 
              title = TeX("$p = 500, \\Delta = 0.5$"))

p2 = plot_acc(folder = "simulation3/patternMOB1", dim = 1000, 
              title = TeX("$p = 1000, \\Delta = 0.5$"))

p3 = plot_acc(folder = "simulation3/patternMOB2", dim = 500, 
              title = TeX("$p = 500, \\Delta = 1$"))

p4 = plot_acc(folder = "simulation3/patternMOB2", dim = 1000, 
              title = TeX("$p = 1000, \\Delta = 1$"))

p5 = plot_acc(folder = "simulation3/patternMOB3", dim = 500, 
              title = TeX("$p = 500, \\Delta = 1.5$"))

p6 = plot_acc(folder = "simulation3/patternMOB3", dim = 1000, legend_pos = c(0.25, 0.15),
              title = TeX("$p = 1000, \\Delta = 1.5$"))

library(cowplot)
p <- plot_grid(p1, p3, p5, p2, p4, p6, byrow = T, nrow = 2, ncol = 3)

ggsave(p, filename = "simulation3/accuracy.pdf", width = 9, height = 6)