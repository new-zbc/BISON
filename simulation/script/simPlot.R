readResult = function(fold, method, dim, col = 1){
  file1 = paste0(fold, "/P_", dim, "/prop_0/summary/", method, ".txt")
  x1 = read.table(file1, sep =" ", header = T)
  x1 = c(mean(x1[, col], na.rm = T), sd(x1[, col], na.rm = T))
  
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
  
  result = rbind(x1, x2, x3, x4, x5)
  colnames(result) = c("ARI", "SD")
  return(result)
}


plot_sim <-function(folder = "simulation1/patternMOB1", dim = 500, col = 1, legend_pos = "none",
                    title = TeX("$p = 500, \\Delta = 1$")){
  method1 = readResult(folder, "bison", dim, col)
  method2 = readResult(folder, "BC", dim, col)
  method3 = readResult(folder, "sparseBC", dim, col)
  method4 = readResult(folder, "spartaco", dim, col)
  method5 = readResult(folder, "Kmeans", dim, col)
  
  df = rbind(method1, method2, method3, method4, method5)
  method = factor(rep(c("Bison", "BC", "sparseBC", "spartaco", "Kmeans"), each = 5), 
                  levels = c("Bison", "BC", "sparseBC", "spartaco", "Kmeans"))
  df = data.frame(df, ratio=rep((0.2*0:4), 5), method = method)
  
  library(ggplot2)
  library(latex2exp)
  
  plotARI <- ggplot(data=df) + 
    geom_line( mapping = aes(x=ratio, y = ARI, linetype=method, color=method)) + 
    geom_point( mapping = aes(x=ratio, y = ARI, color=method, shape=method)) +
    ylim(-0.05, 1.05) + 
    scale_linetype_manual(values = c(5, 4, 3, 2, 1)) +
    scale_color_manual(values = c('purple','darkblue',"orange", "darkred", "grey")) +
    scale_shape_manual(values = c(21, 22, 23, 24, 25, 26)) + 
    geom_errorbar(aes(ymin = ARI-SD, ymax = ARI + SD,x=ratio, y = ARI, color=method), width = 0.01)+
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

p1 = plot_sim(folder = "simulation3/patternMOB1", dim = 500, col = 1,
              title = TeX("$p = 500, \\Delta = 0.5$"))

p2 = plot_sim(folder = "simulation3/patternMOB1", dim = 1000, col = 1,
              title = TeX("$p = 1000, \\Delta = 0.5$"))

p3 = plot_sim(folder = "simulation3/patternMOB2", dim = 500, col = 1,
              title = TeX("$p = 500, \\Delta = 1$"))

p4 = plot_sim(folder = "simulation3/patternMOB2", dim = 1000, col = 1,
              title = TeX("$p = 1000, \\Delta = 1$"))

p5 = plot_sim(folder = "simulation3/patternMOB3", dim = 500, col = 1,
              title = TeX("$p = 500, \\Delta = 1.5$"))

p6 = plot_sim(folder = "simulation3/patternMOB3", dim = 1000, col = 1, legend_pos = c(0.25, 0.2),
              title = TeX("$p = 1000, \\Delta = 1.5$"))

library(patchwork)
library(cowplot)
p <- plot_grid(p1, p3, p5, p2, p4, p6, byrow = T, nrow = 2, ncol = 3)

ggsave(p, filename = "simulation3/ARI_col.pdf", width = 9, height = 6)



################################
#clustering result of row
##################################
p1 = plot_sim(folder = "simulation3/patternMOB1", dim = 500, col = 2,
              title = TeX("$p = 500, \\Delta = 0.5$"))

p2 = plot_sim(folder = "simulation3/patternMOB1", dim = 1000, col = 2,
              title = TeX("$p = 1000, \\Delta = 0.5$"))

p3 = plot_sim(folder = "simulation3/patternMOB2", dim = 500, col = 2,
              title = TeX("$p = 500, \\Delta = 1$"))

p4 = plot_sim(folder = "simulation3/patternMOB2", dim = 1000, col = 2,
              title = TeX("$p = 1000, \\Delta = 1$"))

p5 = plot_sim(folder = "simulation3/patternMOB3", dim = 500, col = 2,
              title = TeX("$p = 500, \\Delta = 1.5$"))

p6 = plot_sim(folder = "simulation3/patternMOB3", dim = 1000, col = 2, legend_pos = c(0.25, 0.2),
              title = TeX("$p = 1000, \\Delta = 1.5$"))

library(patchwork)
p <- plot_grid(p1, p3, p5, p2, p4, p6, byrow = T, nrow = 2, ncol = 3)

ggsave(p, filename = "simulation3/ARI_row.pdf", width = 9, height = 6)


