#' beta strand heatmap
#' 
#' This function allows you to plot beta strand heatmap.
#' @param ddG free energy data
#' @param assay_sele assay_sele
#' @return Nothing
#' @export
Plot_beta_strand_propagation_heattmap<-function(ddG=ddG,
                                                assay_sele=assay_sele){
  ddG<-Read_ddG(ddG=ddG,assay_sele = assay_sele)
  aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))
  heatmap_tool<-data.table(wt_codon = rep(unlist(strsplit(wt_aa,"")),each=20),
                           Pos_real = rep(2:188,each=20),
                           mt_codon = unlist(aa_list))
  ddG_merge<-merge(ddG,heatmap_tool,by=c("wt_codon","Pos_real","mt_codon"),all=T)
  plot_df<-ddG_merge[Pos_real%in%c(35:41,54:60,4:10,77:83,111:117,140:143),]
  plot_df[,pos_wtcodon:=paste0(wt_codon,Pos_real)]
  plot_df[wt_codon==mt_codon,`mean_kcal/mol`:=0]
  plot_df<-within(plot_df, 
                  mt_codon <- factor(mt_codon, 
                                     levels = c("D","E","R","H","K","S","T","N","Q","C","G","P","A","V","I","L","M","F","W","Y")))
  
  plot_df<-within(plot_df, 
                  pos_wtcodon <- factor(pos_wtcodon, 
                                        levels = c("R41","Y40","S39","D38","E37","I36","T35",
                                                   "D54","I55","L56","D57","T58","A59","G60",
                                                   "Y4","K5","L6","V7","V8","V9","G10",
                                                   "G77","F78","L79","C80","V81","F82","A83",
                                                   "M111","V112","L113","V114","G115","N116","K117",
                                                   "P140","F141","I142","E143")))
  ggplot(plot_df,
         aes(y=assay,x=mt_codon))+
    geom_tile(aes(fill=`mean_kcal/mol`))+
    geom_text(data=plot_df[wt_codon==mt_codon,],
              mapping=aes(y=assay,x=mt_codon,label="-"),size=5*0.35)+
    scale_fill_gradient2(low=colour_scheme[["blue"]],mid="gray",high=colour_scheme[["red"]],na.value = "white")+
    labs(fill="\u0394\u0394G(kcal/mol)",)+
    xlab("Mutant aa")+
    ylab(NULL)+
    facet_wrap(~pos_wtcodon,ncol=7)+
    theme_classic2()+
    theme(text = element_text(size=7),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          legend.text = element_text(size=5),
          axis.text.x = element_text(family="Courier New",size =5, vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
          legend.key.height= unit(3.1, 'mm'),
          legend.key.width = unit(3.1, 'mm'),
          legend.key.size = unit(1,"mm"),
          panel.spacing.y = unit(1, "mm"),
          panel.spacing.x = unit(1, "mm"),
          plot.margin=margin(0,0,0,0),
          strip.background = element_rect(color="black",fill=NULL,linetype="solid"))+
    coord_fixed()
}