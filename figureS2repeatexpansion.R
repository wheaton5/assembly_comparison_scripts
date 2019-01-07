region_start = 13091592
region_end = 13627685
p9 = ggplot(subset(dotplot,primary_chrom=="2L" & read_name=="000066F_arrow_arrow"))+
  geom_segment(aes(x=ref_pos_start/1000000,xend=ref_pos_end/1000000,
  y=(read_pos_start+primary_chrom_offset)/1000000,
  yend=(read_pos_end+primary_chrom_offset)/1000000),
  size=0.5,lineend="butt")+
  geom_hline(aes(yintercept=primary_chrom_offset/1000000),size=0.2,linetype="dashed")+
  geom_vline(data=subset(PESTy,(primary_chrom=="2L" ) & primary_chrom!="Mt" & primary_chrom!="UNKN"),
             aes(xintercept=stop/1000000),size=0.2,linetype="dashed")+
  #facet_wrap(~primary_chrom)+
  theme_bw(base_size=12)+
  xlab("")+
  ylab("")+
  theme(panel.grid.minor.x = element_blank() ,
        panel.grid.minor.y = element_blank())+#,
  #panel.grid.minor.x = element_blank() ,
  #panel.grid.minor.y = element_blank())+
  #coord_cartesian(xlim=c(region_start/1000000,region_end/1000000),
  #                ylim=c(2.055e7/1000000,2.1125e7/1000000))+
  scale_x_continuous(limits=c(region_start/1000000,region_end/1000000))+
  scale_y_continuous(limits=c(2.055e7/1000000,2.1125e7/1000000))+
  geom_segment(data=subset(n_s_bed,X1=="2L"),aes(x=X2/1000000,xend=X3/1000000,y=2.06e7/1000000,yend=2.06e7/1000000),size=12)+
  annotation_custom(
    grob = textGrob(label = "Ns track",hjust=0),
    ymin = 2.06e7/1000000,      # Vertical position of the textGrob
    ymax = 2.06e7/1000000,
    xmin = 1.3475e7/1000000,         # Note: The grobs are positioned outside the plot area
    xmax = 1.3475e7/1000000)
gt <- ggplot_gtable(ggplot_build(p9))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
#grid.draw(gt)

p8 = ggplot(covrep)+geom_line(aes(x=(X1)/1000000,y=X2),size=0.3)+
  theme_bw(base_size=12)+
  scale_y_continuous(breaks=c(0,500,1000))+
  scale_x_continuous(limits=c(region_start/1000000,region_end/1000000))+
  geom_hline(aes(yintercept=225),size=0.3)+
  xlab("PEST reference position in Mb on chrom 2L")+
  ylab("coverage")
xoffset = 20550000
p7 = ggplot(cov66F)+geom_line(aes(x=(X1-5000+xoffset)/1000000,y=X2),size=0.3)+
  theme_bw(base_size=12)+
  scale_y_continuous(breaks=c(0,250,500))+
  scale_x_continuous(limits=c(0+xoffset/1000000,(570000+xoffset)/1000000),breaks=c(20600000/1000000,20800000/1000000,21000000/1000000))+
  geom_hline(aes(yintercept=225),size=0.3)+
  xlab("contig 66F position plus offset in Mb")+
  ylab("coverage")+
  #theme(axis.title.y = element_text(angle = 0))+
  coord_flip()

grid.arrange(gt,p8,p7, nrow = 2,
             heights = c(5,1.5),widths=c(1.5,.025,0.34,6),
             layout_matrix=rbind(c(3,3,1,1),c(NA,2,2,2)))



