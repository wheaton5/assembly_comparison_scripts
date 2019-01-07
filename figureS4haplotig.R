
boxx = 23.5
boxxend = 26
boxy = -.5
boxyend = 2.6
basesize=8
boxoutline = 0.4
p13 = ggplot()+
  geom_segment(data = subset(coverage_curated_regions,X1!="Y_unplaced"),
               aes(x=X2/1000000,xend=X3/1000000,y=X4,yend=X4),size=2.5)+
  geom_segment(data = subset(n_s_bed,X1!="Y_unplaced"),
               aes(x=X2/1000000,xend=X3/1000000,y=0,yend=0),size=2.5,color="red")+
  theme_bw(base_size=basesize)+
  scale_y_continuous(limits=c(-.5,4))+
  facet_wrap(~chrom,ncol=1,strip.position="right")+
  xlab("PEST reference position in Mb")+
  ylab("assembly contig alignments to PEST coverage")+
  geom_segment(data = data.frame(chrom="3L",x=boxx,xend=boxxend,y=boxy,yend=boxy),
               aes(x=x,xend=xend,y=y,yend=yend),size=boxoutline)+
  geom_segment(data = data.frame(chrom="3L",x=boxx,xend=boxxend,y=boxyend,yend=boxyend),
               aes(x=x,xend=xend,y=y,yend=yend),size=boxoutline)+
  geom_segment(data = data.frame(chrom="3L",x=boxx,xend=boxx,y=boxy,yend=boxyend),
               aes(x=x,xend=xend,y=y,yend=yend),size=boxoutline)+
  geom_segment(data = data.frame(chrom="3L",x=boxxend,xend=boxxend,y=boxy,yend=boxyend),
               aes(x=x,xend=xend,y=y,yend=yend),size=boxoutline)


refstart = 23500000
refend = 26000000
chrom_of_interest = "3L"
p14 = ggplot(subset(dotplot,
                    primary_chrom==chrom_of_interest &
                      (read_name=="000057F_arrow_arrow" | 
                         read_name == "000050F_arrow_arrow" | 
                         read_name == "000021F_arrow_arrow")))+

  geom_segment(aes(x=ref_pos_start/1000000,xend=ref_pos_end/1000000,
    y=(read_pos_start+primary_chrom_offset)/1000000,
    yend=(read_pos_end+primary_chrom_offset)/1000000, color=read_name),
    size=0.5,lineend="butt")+
  geom_hline(aes(yintercept=primary_chrom_offset/1000000),size=0.2,linetype="dashed")+
  geom_vline(data=subset(PESTy,(primary_chrom==chrom ) & primary_chrom!="Mt" & primary_chrom!="UNKN"),
             aes(xintercept=stop/1000000),size=0.2,linetype="dashed")+
  #facet_wrap(~primary_chrom)+
  theme_bw(base_size=basesize)+
  xlab("PEST reference position in Mb on chrom 3L")+
  ylab("contig position plus offset in Mb")+
  theme(panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank() ,
        panel.grid.minor.y = element_blank())+
  #geom_segment(aes(x=23.0,xend=24068000/1000000,y=23990000/1000000,yend=23990000/1000000),size=0.2)+
  #geom_hline(aes(yintercept=23990000/1000000),size=0.2)+
  geom_rect(data=data.frame(xmin=23.0,xmax=24068000/1000000,ymin=23.8,ymax=23990000/1000000),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),size=0,fill="black",alpha=0.2)+
  geom_rect(data=data.frame(xmin=23.0,xmax=24864000/1000000+.35,ymin=24685000/1000000,ymax=25.035),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),size=0,fill="black",alpha=0.2)+
  geom_rect(data=data.frame(xmin=24068000/1000000-.206,xmax=24068000/1000000,ymin=23,ymax=23.8),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),size=0,fill="black",alpha=0.2)+
  geom_rect(data=data.frame(xmin=24864000/1000000,xmax=24864000/1000000+.35,ymin=23,ymax=24685000/1000000),aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),size=0,fill="black",alpha=0.2)+
  
  #geom_segment(aes(x=24068000/1000000-.206,xend=24068000/1000000-.206,y=23,yend=23.8),size=0.18)+
  #geom_segment(aes(x=23.0,xend=24864000/1000000,y=24685000/1000000,yend=24685000/1000000),size=0.2)+
  #geom_hline(aes(yintercept=24685000/1000000),size=0.2)+
  #geom_segment(aes(x=24068000/1000000,xend=24068000/1000000,y=23,yend=23990000/1000000),size=0.2)+
  #geom_segment(aes(x=24864000/1000000+.35,xend=24864000/1000000+.35,y=23,yend=25.035),size=.2)+
  #geom_vline(aes(xintercept=24068000/1000000),size=0.2)+
  #geom_segment(aes(x=24864000/1000000,xend=24864000/1000000,y=23,yend=24685000/1000000),size=0.2)+
  #geom_vline(aes(xintercept=24864000/1000000),size=0.2)+
  scale_color_manual(breaks=
    c("000021F_arrow_arrow",
      "000057F_arrow_arrow",
      "000050F_arrow_arrow"),
    labels=c("21F","57F","50F"),values=c("red","green3","blue"))+
  guides(color=guide_legend(title="contig")) +
  #scale_x_continuous(limits=c(refstart/1000000,refend/1000000))+
  coord_cartesian(xlim=c(refstart/1000000,refend/1000000),ylim=c(refstart/1000000,refend/1000000))#+
  #scale_y_continuous(limits=c(refstart/1000000,refend/1000000))#+
  #geom_segment(data=subset(n_s_bed,X1==chrom),aes(x=X2,xend=X3,y=1.8e7,yend=1.8e7),size=12)+
  #annotation_custom(
  #  grob = textGrob(label = "Ns track",hjust=0),
  #  ymin = 1.8e7,      # Vertical position of the textGrob
  #  ymax = 1.8e7,
  #  xmin = 2.3e7,         # Note: The grobs are positioned outside the plot area
  #  xmax = 2.3e7)
#gt <- ggplot_gtable(ggplot_build(p14))
#gt$layout$clip[gt$layout$name == "panel"] <- "off"
#grid.draw(gt)
xoffset = 23810000
p15 = ggplot(cov57F)+
  geom_line(aes(x=(X1-5000+xoffset)/1000000,y=X2),size=0.3,color="blue")+
  theme_bw(base_size=basesize)+
  scale_y_continuous(breaks=c(0,150,300))+
  scale_x_continuous()+
  geom_hline(aes(yintercept=225),size=0.2)+
  geom_rect(data=data.frame(xmin=23,xmax=23990000/1000000,ymin=-50,ymax=650),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            size=0,fill="black",alpha=0.2)+
  geom_rect(data=data.frame(xmin=24685000/1000000,xmax=26,ymin=-50,ymax=650),
            aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            size=0,fill="black",alpha=0.2)+
  coord_cartesian(xlim=c(23860000/1000000,24950000/1000000),ylim=c(0,450))+
  xlab("contig 57F position plus offset in Mb")+
  #geom_vline(aes(xintercept=23990000/1000000),size=0.2)+
  #geom_vline(aes(xintercept=24685000/1000000),size=0.2)+
  ylab("coverage")

grid.arrange(p13, p14,p15, nrow = 3,heights = c(5,5,1.7),widths=c(1.56,5,3.91),
             layout_matrix=rbind(c(1,1,1),c(2,2,2),c(NA,3,NA)))

