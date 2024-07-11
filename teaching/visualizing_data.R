# Visualizing your data
# A workshop
# Written by Chelsea Wood (chelwood@uw.edu)
# 11 July 2024


# What is the difference between a continuous and a categorical variable?


# Let us begin by loading the dataset we want to work with

pim_vig_data<-read.csv("data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv")


# A plot with continuous x and continuous y (Did psite_count change over time?)

pim_vig_plot_1<-ggplot(pim_vig_data,aes(YearCollected,psite_count))+
  #scale_color_manual(values=c("#0571b0","#ca0020"))+
  geom_point(size=4)+
  #geom_errorbar(data=pim_vig_data,mapping=aes(x=x,ymin=conf.low,ymax=conf.high),width=0.03)+
  #geom_line(aes(group=group))+
  xlab("Year collected")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
#scale_x_discrete(limits=rev(levels(pr_port_dt$x)),labels=c("historical","contemporary"))+
#theme(legend.position="none")


# A plot with continuous x and continuous y and multiple series (Was change in psite_count over time
# different in control versus impact sites?)

pim_vig_plot_2<-ggplot(pim_vig_data,aes(YearCollected,psite_count,color=CI),grouping=CI,color=CI)+
  #scale_color_manual(values=c("#0571b0","#ca0020"))+
  geom_point(size=4)+
  #geom_errorbar(data=pim_vig_data,mapping=aes(x=x,ymin=conf.low,ymax=conf.high),width=0.03)+
  #geom_line(aes(group=group))+
  xlab("Latitude (degrees N)")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
#scale_x_discrete(limits=rev(levels(pr_port_dt$x)),labels=c("historical","contemporary"))+
#theme(legend.position="none")


# A plot with categorical x and continuous y (Was psite_count different between control and impact sites?)

pim_vig_plot_3<-ggplot(pim_vig_data,aes(CI,psite_count))+
  #scale_color_manual(values=c("#0571b0","#ca0020"))+
  geom_boxplot(size=4)+
  #geom_errorbar(data=pim_vig_data,mapping=aes(x=x,ymin=conf.low,ymax=conf.high),width=0.03)+
  #geom_line(aes(group=group))+
  xlab("")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
  #scale_x_discrete(limits=rev(levels(pr_port_dt$x)),labels=c("historical","contemporary"))+
  #theme(legend.position="none")


# A plot with categorical x and continuous y and multiple series (Was psite_count different between control and
# impact sites, and was this difference consistent across decades?)

pim_vig_plot_4<-ggplot(pim_vig_data,aes(CI,psite_count,color=combo),grouping=CI,color=combo)+
  #scale_color_manual(values=c("#0571b0","#ca0020"))+
  geom_boxplot(size=4)+
  #geom_errorbar(data=pim_vig_data,mapping=aes(x=x,ymin=conf.low,ymax=conf.high),width=0.03)+
  #geom_line(aes(group=group))+
  xlab("")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
#scale_x_discrete(limits=rev(levels(pr_port_dt$x)),labels=c("historical","contemporary"))+
#theme(legend.position="none")
