#' all ddG point plot
#' 
#' This function allows you to plot all ddG heatmap.
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
#' @param ddG7 free energy data7
#' @param assay7 data7's name
#' @param anno annotation df final_distance_dc_anno_for_anno_mt_4
#' @param rect_input rect for beta strands
#' @param rect_alpha rect for alpha 1
#' @return Nothing
#' @export
Plot_figure4A_point_all_ddG_alphabeta<-function(
  ddG1=ddG1,
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
  ddG7=ddG7,
  assay7=assay7,
  anno=anno,
  rect_input=rect_input,
  rect_alpha=rect_alpha
){
  ddG1<-Read_ddG(ddG1,assay1)
  ddG2<-Read_ddG(ddG2,assay2)
  ddG3<-Read_ddG(ddG3,assay3)
  ddG4<-Read_ddG(ddG4,assay4)
  ddG5<-Read_ddG(ddG5,assay5)
  ddG6<-Read_ddG(ddG6,assay6)
  ddG7<-Read_ddG(ddG7,assay7)
  rects_dt<-as.data.table(rect_input)
  rects_alphas<-as.data.table(rect_alpha)
  all_ddG<-rbind(ddG1,ddG2,ddG3,ddG4,ddG5,ddG6,ddG7)
  all_data_25_anno_bs<-merge(all_ddG,
                             anno,
                             by.x=c("Pos_real","assay"),
                             by.y=c("Pos","assay"),all=T)
  all_data_25_anno_bs[is.na(binding),binding:="no"]
  all_data_25_anno_bs<-within(all_data_25_anno_bs, 
                              assay <- factor(assay,
                                              levels=c("folding",
                                                       "RAF1",
                                                       "PIK3CG",
                                                       "RALGDS",
                                                       "SOS1",
                                                       "DARPin K27",
                                                       "DARPin K55")))
  all_data_25_anno_bs<-within(all_data_25_anno_bs, 
                              binding <- factor(binding,
                                                levels=c("core",
                                                         "surface",
                                                         "HVR",
                                                         "GTP",
                                                         "GDP",
                                                         "both",
                                                         "RAF",
                                                         "PI3",
                                                         "RAL",
                                                         "SOS",
                                                         "K27",
                                                         "K55",
                                                         "no")))
  ggplot()+
    geom_rect(data=rects_dt,aes(ymin=-2,ymax=3.5,xmin=xstart-0.5,xmax=xend+0.5),fill="black",alpha=0.08)+
    geom_rect(data=rects_alphas,aes(ymin=-2,ymax=3.5,xmin=xstart-0.5,xmax=xend+0.5),fill="black",alpha=0.05)+
    geom_jitter(data=all_data_25_anno_bs[id!="WT"],aes(x=Pos_real,y=`mean_kcal/mol`,color=binding),size=0.1,width=0.4,height = 0)+
    facet_wrap(~assay,nrow = 7)+
    scale_x_continuous(expand=c(1/188,11/188),breaks=seq(0,187,5),
                       minor_breaks = seq(2,187,1))+
    scale_color_manual(values=c(colour_scheme[["dark green"]], colour_scheme[["orange"]],"gray",
                                colour_scheme[["blue"]],colour_scheme[["blue"]],
                                colour_scheme[["purple"]],colour_scheme[["red"]],colour_scheme[["red"]],
                                colour_scheme[["red"]],colour_scheme[["red"]],colour_scheme[["red"]],colour_scheme[["red"]],"gray"))+
    theme_classic()+
    theme(text = element_text(size = 7),
          axis.text = element_text(size = 7),legend.text = element_text(size = 7),
          legend.position="bottom",
          panel.spacing.y = unit(1, "mm"),
          strip.background = element_rect(color="black",fill=NULL,linetype="solid")
    )+
    coord_fixed(ratio=2)
}