#' Plot major allosteric sites
#' 
#' This function allows you to plot major allosteric sites.
#' @param ddG1 free energy data1
#' @param assay1 data1's name
#' @param ddG2 free energy data2
#' @param assay2 data2's name
#' @param ddG3 free energy data3
#' @param assay3 data3's name
#' @param ddG4 free energy data4
#' @param assay4 data4's name
#' @param ddG5 free energy data5
#' @param assay5 data5's name
#' @param ddG6 free energy data6
#' @param assay6 data6's name
#' @param anno anno_final3
#' @param allosteric_list allosteric_list
#' @param rect_input rect_input
#' @return Nothing
#' @export
#' 
Plot_major_allosteric_site_combined_threshold<-function(ddG1=ddG1,
                                                             assay1=assay1,
                                                             ddG2=ddG2,
                                                             assay2=assay2,
                                                             ddG3=ddG3,
                                                             assay3=assay3,
                                                             ddG4=ddG4,
                                                             assay4=assay4,
                                                             ddG5=ddG5,
                                                             assay5=assay5,
                                                             ddG6=ddG6,
                                                             assay6=assay6,
                                                             anno=anno_final3,
                                                             allosteric_list = allosteric_list,
                                                             rect_input=rect_input){
  weighted_mean_ddG<-list()
  weighted_mean_ddG[[assay1]]<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG1,
                                                              assay_sele = assay1)
  weighted_mean_ddG[[assay2]]<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG2,
                                                              assay_sele = assay2)
  weighted_mean_ddG[[assay3]]<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG3,
                                                              assay_sele = assay3)
  weighted_mean_ddG[[assay4]]<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG4,
                                                              assay_sele = assay4)
  weighted_mean_ddG[[assay5]]<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG5,
                                                              assay_sele = assay5)
  weighted_mean_ddG[[assay6]]<-Get_weighted_mean_abs_ddG_mutcount(ddG=ddG6,
                                                              assay_sele = assay6)
  
  assay_list<-list(assay1,assay2,assay3,assay4,assay5,assay6)
  
  data_plot<-data.table()
  for(assayi in assay_list){
    data_plot_assayi<-merge(weighted_mean_ddG[[assayi]],anno,by="Pos",all=T)
    data_plot_assayi[,binding_type:="allosteric site"]
    data_plot_assayi[get(paste0("scHAmin_ligand_",assayi))<5,binding_type:="binding site"]
    data_plot_assayi[,binding_type_gtp_included:=binding_type]
    data_plot_assayi[get(paste0("GXPMG_scHAmin_ligand_",assayi))<5,binding_type_gtp_included:="GTP binding site"]
    data_plot<-rbind(data_plot,data_plot_assayi)
  }
  reg_threshold<-data_plot[binding_type=="binding site",sum(abs(.SD[[1]])/.SD[[2]]^2, na.rm = T)/sum(1/.SD[[2]]^2, na.rm = T),.SDcols = c("mean","sigma")]
  
  
  data_plot[,site_type:="Reminder"]
  data_plot[binding_type_gtp_included=="binding site",site_type:="Binding interface site"]
  data_plot[binding_type_gtp_included=="GTP binding site",site_type:="Other GTP pocket site"]
  data_plot[binding_type_gtp_included=="GTP binding site"&mean>reg_threshold&binding_type!="binding site"&count>9.5,site_type:="Allosteric GTP pocket site"]
  data_plot[binding_type_gtp_included=="allosteric site"&mean>reg_threshold&count>9.5,site_type:="Major allosteric site"]
  data_plot[,colors_type:="others"]
  rects_dt<-as.data.table(rect_input)
  data_plot[Pos>=rects_dt[col=="b1",xstart]&Pos<=rects_dt[col=="b1",xend],colors_type:="b1"]
  data_plot[Pos>=rects_dt[col=="b2",xstart]&Pos<=rects_dt[col=="b2",xend],colors_type:="b2"]
  data_plot[Pos>=rects_dt[col=="b3",xstart]&Pos<=rects_dt[col=="b3",xend],colors_type:="b3"]
  data_plot[Pos>=rects_dt[col=="b4",xstart]&Pos<=rects_dt[col=="b4",xend],colors_type:="b4"]
  data_plot[Pos>=rects_dt[col=="b5",xstart]&Pos<=rects_dt[col=="b5",xend],colors_type:="b5"]
  data_plot[Pos>=rects_dt[col=="b6",xstart]&Pos<=rects_dt[col=="b6",xend],colors_type:="b6"]
  data_plot[,shape_beta:="others"]
  data_plot[colors_type!="others",shape_beta:="beta strand"]
  data_plot<-data_plot[Pos_real>1&count>9.5,]
  data_plot<-within(data_plot,
                    site_type<-factor(site_type,
                                      levels = c("Binding interface site",
                                                 "Allosteric GTP pocket site",
                                                 "Other GTP pocket site",
                                                 "Major allosteric site",
                                                 "Reminder")))
  data_plot<-within(data_plot,
                    assay<-factor(assay,
                                  levels = c("RAF",
                                             "PI3",
                                             "RAL",
                                             "SOS",
                                             "K27",
                                             "K55")))
  allosteric_list<-list()
  for(assayi in assay_list){
    data_plot[assay==assayi,distance_bp:=get(paste0("scHAmin_ligand_",assayi))]
    allosteric_list[[assayi]][["Binding interface site"]]<-data_plot[binding_type=="binding site"&assay==assayi,Pos]
    allosteric_list[[assayi]][["Allosteric GTP pocket site"]]<-data_plot[site_type=="Allosteric GTP pocket site"&assay==assayi,Pos]
    allosteric_list[[assayi]][["Other GTP pocket site"]]<-data_plot[site_type=="Other GTP pocket site"&assay==assayi,Pos]
    allosteric_list[[assayi]][["Major allosteric site"]]<-data_plot[site_type=="Major allosteric site"&assay==assayi,Pos]
  }
  
  allosteric_list[["p"]]<-ggplot()+
    geom_point(data=data_plot,
               mapping = aes(x=distance_bp,
                             y=mean,
                             color=site_type,
                             shape=as.factor(shape_beta)),
               size=0.35)+
    scale_color_manual(values=c(colour_scheme[["red"]],
                                colour_scheme[["blue"]],
                                colour_scheme[["light blue"]],
                                colour_scheme[["green"]],
                                "gray"),
                       labels=c("Binding interface",
                                "Allosteric GTP pocket",
                                "Other GTP pocket",
                                "Major allosteric",
                                "Others"
                       ))+
    geom_pointrange(data=data_plot,aes(x=distance_bp,
                                       y=mean,
                                       color=site_type,ymin=mean-sigma, ymax=mean+sigma,
                                       shape=as.factor(shape_beta)),size=0.35)+
    geom_hline(yintercept = reg_threshold,linetype =2,size=0.1)+
    geom_vline(xintercept = 5,linetype =2,size=0.1)+
    geom_hline(yintercept = 0,linetype ="solid",size=0.1)+
    geom_vline(xintercept = 0,linetype ="solid",size=0.1)+
    geom_text_repel(data=data_plot[site_type=="Major allosteric site",],aes(x=distance_bp,
                                                                            y=mean,label=Pos),nudge_y=0.05,color=colour_scheme[["green"]],size=5*0.35)+
    geom_text_repel(data=data_plot[site_type=="Allosteric GTP pocket site",],aes(x=distance_bp,
                                                                                 y=mean,label=Pos),nudge_y=0.05,color=colour_scheme[["blue"]],size=5*0.35)+
    xlab(expression(paste("Distance to binding parnter ("*ring(A)*")")))+
    ylab("Weighted mean |ddG| (kcal/mol)")+
    labs(color=NULL)+
    facet_wrap(~assay,ncol=3)+
    theme_classic()+
    theme(axis.text.x = element_text(size =7, vjust=.5, hjust=.5),
          axis.text.y = element_text(size=7, vjust = .5,hjust =.5),
          text = element_text(size=7),
          legend.position="right",
          strip.text.x = element_text(size=7),
          strip.background = element_rect(colour="white", fill="white"),
          panel.spacing = unit(0.2, "mm"),
          #legend.key.height= unit(3.1, 'mm'),
          #legend.key.width = unit(5, 'mm'),
          legend.text = element_text(size=5),plot.margin=margin(0,1,0,1,"mm"),
          legend.margin=margin(0,0,0,l=-2,"mm"),
          legend.spacing.y = unit(0, 'mm'),
          legend.key.height=unit(4,"mm"))+
    labs(shape=NULL)
  return(allosteric_list)
}