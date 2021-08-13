################################################################################
### combined interference plot script
################################################################################

### load packages

library(ggplot2); library(cowplot); library(RColorBrewer);

### load data for plot

load(file = 'Female_Int_Plot.RData')

load(file = 'Male_Int_Plot.RData')

load(file = 'Juvenile_Int_Plot.RData')

### make plots for each sex/stage

Female <- ggplot(data = Female_Int_DataFrame, aes(x = Density, y = Q50/24, color = SexStage)) + geom_line(size = 1) + 
  geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24, fill = SexStage), alpha = 0.25, color = NA) + 
  geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24, fill = SexStage), alpha = 0.5, color = NA) + 
  theme_cowplot() + xlab('Density of Sex/Stage') + ylab('Feeding Rate \n(Midges per hour)') + labs(fill = 'Sex/Stage') +
  labs(color = 'Sex/Stage') + scale_color_brewer(palette = 'Dark2') + scale_fill_brewer(palette = 'Dark2') + ggtitle('Females') 

Male <- ggplot(data = Male_Int_DataFrame, aes(x = Density, y = Q50/24, color = SexStage)) + geom_line(size = 1) + 
  geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24, fill = SexStage), alpha = 0.25, color = NA) + 
  geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24, fill = SexStage), alpha = 0.5, color = NA) + 
  theme_cowplot() + xlab('Density of Sex/Stage') + ylab('Feeding Rate \n(Midges per hour)') + labs(fill = 'Sex/Stage') +
  labs(color = 'Sex/Stage') + scale_color_brewer(palette = 'Dark2') + scale_fill_brewer(palette = 'Dark2') + ggtitle('Males') 

Juvenile <- ggplot(data = Juvenile_Int_DataFrame, aes(x = Density, y = Q50/24, color = SexStage)) + geom_line(size = 1) + 
  geom_ribbon(aes(ymin = Q5/24, ymax = Q95/24, fill = SexStage), alpha = 0.25, color = NA) + 
  geom_ribbon(aes(ymin = Q25/24, ymax = Q75/24, fill = SexStage), alpha = 0.5, color = NA) + 
  theme_cowplot() + xlab('Density of Sex/Stage') + ylab('Feeding Rate \n(Midges per hour)') + labs(fill = 'Sex/Stage') +
  labs(color = 'Sex/Stage') + scale_color_brewer(palette = 'Dark2') + scale_fill_brewer(palette = 'Dark2') + ggtitle('Juveniles')

### combine plots into single plot

Combined <- plot_grid(Female, Male, Juvenile, nrow = 1, ncol = 3)

#save_plot(filename = 'CombinedInterferencePlot.png', plot = Combined, nrow =1, ncol = 3, base_asp = 1.3)

