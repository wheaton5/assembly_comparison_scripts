

p = ggplot(subset(dotplot,primary_chrom=="X" & 
                    (read_name == "000022F_arrow_arrow" |
                       read_name == "000070F_arrow_arrow" |
                       read_name == "000052F_arrow_arrow" | 
                       read_name == "000266F_arrow_arrow")))+
  geom_segment(aes(x=ref_pos_start/1000000,xend=ref_pos_end/1000000,
    y=(read_pos_start+primary_chrom_offset)/1000000,
    yend=(read_pos_end+primary_chrom_offset)/1000000, 
          color=read_name),
    size=0.5,lineend="butt")+
  geom_hline(aes(yintercept=primary_chrom_offset/1000000),size=0.2,linetype="dashed")+
  geom_vline(data=subset(PESTy,(primary_chrom=="X" ) & primary_chrom!="Mt" & primary_chrom!="UNKN"),
             aes(xintercept=stop/1000000),size=0.1,linetype="dashed")+
  #facet_wrap(~primary_chrom)+
  theme_bw(base_size=14)+
  xlab("PEST reference position in MB in chrom X")+
  ylab("")+
  scale_color_discrete(breaks=
                         c("000052F_arrow_arrow",
                           "000266F_arrow_arrow",
                           "000022F_arrow_arrow",
                           "000070F_arrow_arrow"),
                       labels=c("52F","266F","22F","70F"))+
  scale_x_continuous(limits=c(1.75e7/1000000,2.25e7/1000000))+
  scale_y_continuous(limits=c(1.8e7/1000000,2.4e7/1000000))+
 guides(color=guide_legend(title="contig")) +
  geom_segment(data=subset(n_s_bed,X1=="X"),aes(x=X2/1000000,xend=X3/1000000,y=18,yend=18),size=12)+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank() ,
        panel.grid.minor.y = element_blank())+
  annotation_custom(
    grob = textGrob(label = "Ns track",hjust=0),
    ymin = 1.8e7/1000000,      # Vertical position of the textGrob
    ymax = 1.8e7/1000000,
    xmin = 2.3e7/1000000,         # Note: The grobs are positioned outside the plot area
    xmax = 2.3e7/1000000)
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

p2 = ggplot(subset(totalcov2,X1+18019231<2.4e7))+
  geom_line(aes(x=(X1+18019231)/1000000,y=X2),size=0.2)+
  geom_hline(aes(yintercept=225),size=0.2)+
  theme_bw(base_size=14)+
  #theme(axis.title.y = element_text(angle = 0))+
  xlab("contig position plus offset in MB")+
  ylab("coverage")+

  scale_y_continuous(limits=c(0,400),breaks=c(0,150,300))+
  coord_flip()



grid.arrange(p2,gt, nrow = 1,widths = c(1,5))