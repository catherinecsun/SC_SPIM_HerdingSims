coverages<-rbind.fill(coverages_SC,coverages_SPIM)
  
ggplot(coverages[coverages$p0==0.05&coverages], aes(x=aggregation, y=coverage, group=cohesion)) +
  geom_line(aes(color=cohesion))+
  geom_point(aes(color=cohesion))+
  facet_grid(rows=vars(p0),cols=vars(param),labeller = labeller(p0 = p0.labs))+
  labs(title="Abundance (N)       and Sigma (\u03c3)", x="Aggregation (Group Size)",y="Coverage")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))

