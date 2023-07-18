###beginning
library(ggpubr)
library(GGally)
#library(bio3d)
library(ggplot2)
library(plotly)
library(ComplexHeatmap)
library(circlize)
library(seqinr)
library(caTools)
library(pROC)
library(readr)
library(knitr)
library(PRROC)
library(rgl)
library(ggVennDiagram)
library(ggstatsplot)
library(plot3D)
library("devtools")
library(roxygen2)
#library(data.table)
library(ROCR)
library(openxlsx)
library(ggrepel)
#library(Cairo)
######
setwd("/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/")
#
library("devtools")
library(roxygen2)
#create("krasddpcams")
setwd("/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/krasddpcams")
detach(package:krasddpcams, unload = TRUE)
document()
setwd("..")
library("krasddpcams")

######
#basic setting
######
options(ggrepel.max.overlaps = Inf)
colour_scheme<- list(
  "blue"="#1B38A6",#rgb(27, 56, 166)
  "red"="#F4270C",#rgb(244, 39, 12)
  "orange"="#F4AD0C",#rgb(244, 173, 12)
  "green"="#09B636",#rgb(9, 182, 54)
  "yellow"="#F1DD10",#rgb(241, 221, 16)
  "purple"="#6D17A0",#rgb(109, 23, 160)
  "pink"="#FFB0A5",#rgb(255, 176, 165)
  "light orange"="#FFE4A5",#rgb(255, 228, 165)
  "light blue"="#9DACE3",#rgb(157, 172, 227)
  "light green"="#97E9AD",#rgb(151, 233, 173)
  "light red"="#FF6A56",#rgb(255, 106, 86)
  "dark red"="#A31300",#rgb(163, 19, 0)
  "dark blue"="#0C226F",#rgb(12, 34, 111)
  "dark green"="#007A20"#rgb(0, 122, 32)
)
aa_list <- as.list(unlist(strsplit("GAVLMIFYWKRHDESTCNQP", "")))

######
#KRAS information
######
wt_aa<-"TEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKDGKKKKKKSKTKCVIM"
###beta strands
rects<-data.frame(xstart=c(2,37,49,77,111,141),
                  xend=c(10,46,58,83,116,143), 
                  col=c("b1","b2","b3","b4","b5","b6"))
###alpha helix 1
rects_alpha<-data.frame(xstart=c(16),
                        xend=c(25
                        ), 
                        col=c("a1"))
###all alpha helix
rects_alpha_all<-data.frame(xstart=c(16,62,87,127, 152),
                            xend=c(25,74,104, 137,166 ), 
                            col=c("a1","a2","a3","a4","a5"))
##
all_distance_anno<-fread("./DATA/all_distance_anno.csv")
all_long_distance<-fread("./DATA/all_long_distance.csv")
all_long_distance_GXPMG<-fread("./DATA/all_long_distance_GXPMG.csv")
final_distance_dc_anno_for_anno_mt_5<-fread("./DATA/final_distance_dc_anno_for_anno_mt_5.csv")
anno_final3<-fread("./DATA/anno_final3.csv")
######
#loading data
######
stability_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_abundance_1_fitness_replicates_fullseq.RData",
                                                 block2_dimsum_df = "./DATA/CW_RAS_abundance_2_fitness_replicates_fullseq.RData",
                                                 block3_dimsum_df = "./DATA/CW_RAS_abundance_3_fitness_replicates_fullseq.RData")
RAF_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")
PI3_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_PI3_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_PI3_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_PI3_3_fitness_replicates_fullseq.RData")
RAL_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_RAL_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_RAL_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_RAL_3_fitness_replicates_fullseq.RData")
SOS_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_SOS_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_SOS_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_SOS_3_fitness_replicates_fullseq.RData")
K27_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_K27_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_K27_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_K27_3_fitness_replicates_fullseq.RData")
K55_nor_df <- Normalize_growthrate_fitness(block1_dimsum_df = "./DATA/CW_RAS_binding_K55_1_fitness_replicates_fullseq.RData",
                                           block2_dimsum_df = "./DATA/CW_RAS_binding_K55_2_fitness_replicates_fullseq.RData",
                                           block3_dimsum_df = "./DATA/CW_RAS_binding_K55_3_fitness_replicates_fullseq.RData")
stab <- stability_nor_df
RAF <- RAF_nor_df
K27 <- K27_nor_df
PI3 <- PI3_nor_df
RAL <- RAL_nor_df
SOS <- SOS_nor_df
K55 <- K55_nor_df
all_data <- Merge_dimsum_df(stab,RAF,K27,PI3,RAL,SOS,K55)
all_data_pos<-Pos_id(all_data,wt_aa)
RAF_nor_df_pos <- Pos_id(RAF_nor_df,wt_aa)

#########################
#### revision fig mut vs wt block1 RAF1RBD data
#########################
#a doubles vs singles fitness
RAF_nor_df<- Normalize_growthrate_fitness_block1(block1_dimsum_df = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData")
RAF_nor_df_pos<-Pos_id(RAF_nor_df,wt_aa)
Get_background_order_fitness(RAF_nor_df_pos,50)
ggarrange(
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"G12C"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"G12D"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"G12V"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"L6H"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"L19P"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"S17N"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"T2K"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"T35S"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"V14S"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"Y40A"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"D38C"),
Plot_fitness_correlation_doublevssingle(RAF_nor_df_pos,"E37G"))
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figuresS2a_fitness_cor_2vs1.pdf", device = cairo_pdf,height = 150,width=273,units = "mm")
#figure legands
#Scatter plot comparing binding fitness of mutation in mutation background or wild-type background. 
#Substitutions in the binding interface are indicated in red.
ggarrange(
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"G12C"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"G12D"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"G12V"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"L6H"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"L19P"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"S17N"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"T2K"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"T35S"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"V14S"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"Y40A"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"D38C"),
  Plot_fitness_correlation_doublevssingle(GAP_nor_df_pos,"E37G"))
Get_background_order_fitness(GAP_nor_df_pos,100)

RAF_nor_df_3blocks<-Normalize_growthrate_fitness(block1_dimsum_df = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/CW_RAS_binding_RAF_1_fitness_replicates_fullseq.RData",
                                                 block2_dimsum_df = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/CW_RAS_binding_RAF_2_fitness_replicates_fullseq.RData",
                                                 block3_dimsum_df = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/CW_RAS_binding_RAF_3_fitness_replicates_fullseq.RData")


RAF_nor_df_3blocks_pos<-Pos_id(RAF_nor_df_3blocks,wt_aa)
#########################
#### revision single mutations tecan test
#########################
#b growthrate correlation between tecan test vs bindingPCA assay
Plot_cor_gr_tecan_binding(tecandata = "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/invitro validation tecan 0615/20230616_tecanA_1.txt",
                          plasmidid = "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/invitro validation tecan 0615/plasmid_plate_0615.csv",
                          single = all_data_pos[Nham_aa<=1,])
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figureS1D_tecantest_fitness_cor.pdf", device = cairo_pdf,height = 40,width=40,units = "mm")
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures1b_tecantest_fitness_cor.pdf", device = cairo_pdf,height = 60,width=60,units = "mm")
# figure legends
# figure S1D Comparison of individually measured growth rates to those fitness from deep sequencing for KRAS
# The red line corresponds to a linear regression model. Pearson’s r is shown. 
#########################
#### revision double mutations tecan test
#########################
Plot_cor_gr_tecan_binding_MoCHi_prediction_double(tecandata = "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/invitro validation tecan 0615/20230616_tecanA_1.txt",
                                                  plasmidid = "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/invitro validation tecan 0615/plasmid_plate_0615.csv",
                                                  ddG1 = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/weights_Folding.txt",
                                                  ddG2 = "/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/weights_Binding_RAF.txt")
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures1c_tecantest_doublemutations_cor.pdf", device = cairo_pdf,height = 60,width=60,units = "mm")
# figure legends
# figure xxx Comparison of individually measured growth rates of double mutations to 
# model-predict binding fraction (linearly correlated to growth rates).
# The red line corresponds to a linear regression model. Pearson’s r is shown. 

#########################
#### revision sup fig2 GAP PCA
#########################
#a heatmap
GAP_nor_df <- Normalize_growthrate_fitness_block1(block1_dimsum_df = "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/GAP_PCA_data/CW_RAS_binding_RAFGAP_1_fitness_replicates_fullseq.RData")
GAP_nor_df_pos<-Pos_id(GAP_nor_df,wt_aa)
Plot_heatmap_fitness_block1(GAP_nor_df_pos)+
  theme(axis.text.y = element_text(family="Courier New",angle=90,size=5, vjust =.5,hjust =.5,margin=margin(0,-0.5,0,0,"mm")))

ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures2a_GAP_block1_fitness.pdf", device = cairo_pdf,height = 120,width=120,units = "mm")
# figure legends
# Heat maps of fitness effects of single aa substitutions for KRAS-RAF1 coexpression GAP from BindingPCA
# White, missing values; - , wild-type aa; X, STOP codon.

#b fitness correlation
Plot_fitness_correlation_block1(GAP_nor_df_pos)
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures2b_GAP_block1_fitness_cor.pdf", device = cairo_pdf,height = 60,width=60,units = "mm")
# figure legends
# Scatter plots showing the reproducibility of binding fitness of block1 estimates from bindingPCA coexpression with GAP. 
# Pearson’s r indicated on the top right corner. 
# Plot_kuriyan_cor<-function(fitness=all_data_pos,kuriyan_list=kuriyan_list,comparedto=comparedto){
#   kuriyan_list<-readRDS(kuriyan_list)
#   kuriyan_long_df_dc<-Get_dc_long_kuriyan_data(kuriyan_list_data=kuriyan_list)
#   kuriyan_all<-Plot_fitness_kuriyan_data(fitness=fitness,
#                                          kuriyan_data=kuriyan_long_df_dc)
#   kuriyan_lm<-lm(nor_fitness~get(comparedto),kuriyan_all)
#   ggplot(kuriyan_all,aes(x=nor_fitness,y=get(comparedto)))+
#     stat_binhex(bins = 50,size=0,color="black") +
#     scale_fill_gradient(low="white",high="black",trans="log10") +
#     geom_smooth(method=lm,se=T,size=0.3,color=colour_scheme[["red"]]) +
#     annotate("text",x=-0.5,y=0.5,,label = paste0("r = ",round(sqrt(summary(kuriyan_lm)$r.squared),3)) ,size=7*0.35)+
#     ylab("Mean E score(Bandaru et al. 2017)")+
#     xlab("Binding fitness")+theme_classic2()+
#     theme(text = element_text(size=7),
#           legend.position="right",
#           legend.text = element_text(size=7),
#           axis.text.x = element_text(size =7,vjust=.5, hjust=.5),
#           axis.text.y = element_text(size=7, vjust = .5,hjust = .5,margin=margin(0,0,0,0,"mm")),
#           legend.key.height= unit(3.1, 'mm'),
#           legend.key.width = unit(3.1, 'mm'),
#           legend.key.size = unit(1,"mm"),
#           plot.margin=margin(0,0,0,0))
# }
# #c fitness correlation with Kuriyan data
# GAP_nor_df_pos_RAF<-GAP_nor_df_pos[,assay:="RAF"]
# Plot_kuriyan_cor(fitness=GAP_nor_df_pos_RAF[block=="block1",],
#                  kuriyan_list="/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/kuriyan_list.RData",
#                  comparedto="KRas173_GAP")+
#   ggtitle("KRas173_GAP block1")
# Plot_kuriyan_cor(fitness=GAP_nor_df_pos_RAF[block=="block1",],
#                  kuriyan_list="/Users/cweng/Desktop/Chenchun_KRAS_paper/Github_codes_updating/DATA/kuriyan_list.RData",
#                  comparedto="HRas188_BaF3")+
#   ggtitle("HRas188_BaF3 block1")
#########################
#### revision sup fig3 full RAF1
#########################
#a fitness heatmap of block1 of fullRAF1 bindingPCA
fullRAF_nor_df <- Normalize_growthrate_fitness_block1(block1_dimsum_df = "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/0628_merge_data/CW_RAS_binding_FULLRAF_1_fitness_replicates.RData")
fullRAF_nor_df_pos<-Pos_id(fullRAF_nor_df,wt_aa)
Plot_heatmap_fitness_block1(fullRAF_nor_df_pos)+
  theme(axis.text.y = element_text(family="Courier New",angle=90,size=5, vjust =.5,hjust =.5,margin=margin(0,-0.5,0,0,"mm")))
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures3a_fullRAF1_block1_fitness.pdf", device = cairo_pdf,height = 120,width=120,units = "mm")
# figure legends
# Heat maps of fitness effects of single aa substitutions for KRAS-full RAF1 from BindingPCA.
# White, missing values; - , wild-type aa; X, STOP codon.

#b fitness correlation across triplicates
Plot_fitness_correlation_block1(fullRAF_nor_df_pos)
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures3b_fullRAF1_block1_fitness_cor.pdf", device = cairo_pdf,height = 60,width=60,units = "mm")
# figure legends
# Scatter plots showing the reproducibility of binding fitness of block1 for KRAS-full RAF1 from BindingPCA. 
# Pearson’s r indicated on the top right corner. 

#c fitness full RAF1 vs RAF1RBD
Plot_fitness_cor_RAF1_full_RBD(fullRAF_nor_df_pos,RAF_nor_df_pos)
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures3c_fullRAF1_RBD_block1_fitness_cor.pdf", device = cairo_pdf,height = 40,width=60,units = "mm")
# figure legends
# Comparisons of binding fitness of single aa substitutions for KRAS-RAF1RBD and for KRAS-full RAF1 from ddPCA.
# Substitutions in RAF1-RBD binding interface are indicated in red, in RAF1-CRD binding interface are indicated in yellow.

#d mochi fitting Andre's codes
Plot_ddG_fullRAF_RBD("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/Mochi/weights_Binding_RAF.txt",
                     "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/Mochi/weights_Binding_FULLRAF.txt")
ggsave("/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fugures//figures3d_fullRAF1_RBD_block1_ddG_cor.pdf", device = cairo_pdf,height = 40,width=60,units = "mm")
# figure legends
# Comparisons of model-inferred binding free energy changes of single aa substitutions for KRAS-RAF1RBD and for KRAS-full RAF1.
# Substitutions in RAF1-RBD binding interface are indicated in red, in RAF1-CRD binding interface are indicated in yellow.
#########################
#### revision sup fig1 cross triplicates validation of MoChi
#########################
report_outpath <- "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/Mochi/norep3_plots"
dir.create(report_outpath)
plot_model_performance_heldoutrep <- function(
  input_dt,
  report_outpath,
  name,
  highlight_colour = "red"
){
  
  #Plot observed versus predicted fitness - not facetted on mutation order - only validation data
  plot_dt <- input_dt[!is.na(fitness_pred),.(observed_fitness = fitness, predicted_fitness = fitness_pred)]
  d <- ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness (Rep. 3)") +
    ggplot2::ylab("Predicted fitness (Reps. 1&2)") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete"), 2), sep="")),], ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()
  d <- d + ggplot2::geom_abline(color = 'red')
  ggplot2::ggsave(file.path(report_outpath, paste0("fitness_observed_predicted_scatter_allmodels_test_", name, ".pdf")), d, width = 4, height = 3, useDingbats=FALSE)
  
}
pred_RAF_orig <- "./DATA/predicted_phenotypes_all.txt"
pred_RAF_norep3 <- "/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/Mochi/predicted_phenotypes_all.txt"
#Load prediction data
pred_dt_orig <- fread(pred_RAF_orig)
pred_dt_norep3 <- fread(pred_RAF_norep3)
#Merge
pred_dt <- merge(
  pred_dt_norep3[,.SD,,.SDcols = names(pred_dt_norep3)[!names(pred_dt_norep3) %in% c("fitness", "sigma")]],
  pred_dt_orig[, .(aa_seq, phenotype, fitness = fitness3_uncorr, sigma = sigma3_uncorr)],
  by = c("aa_seq", "phenotype"))[!is.na(fitness)]
pred_dt[, fitness_pred := mean]
plot_dt <- pred_dt[!is.na(fitness_pred)&phenotype==1,.(observed_fitness = fitness, predicted_fitness = fitness_pred)]
#Plot performance on held out replicate (blocks separately)
for(i in names(pred_dt)[grepl("Abundance|Binding", names(pred_dt))]){
  plot_model_performance_heldoutrep(
    input_dt = pred_dt[get(i)==1],
    report_outpath = report_outpath,
    name = i)
}

pred_dt[,phenotype]
Plot_pre_12_ob_3_fitness_cor<-function(prediction=prediction,
                                       observed=observed,
                                       phenotypen=phenotypen){
  pred_RAF_orig <- fread(observed)
  pred_dt_norep3 <- fread(prediction)
  pred_dt <- merge(
    pred_dt_norep3[,.SD,,.SDcols = names(pred_dt_norep3)[!names(pred_dt_norep3) %in% c("fitness", "sigma")]],
    pred_dt_orig[, .(aa_seq, phenotype, fitness = fitness3_uncorr, sigma = sigma3_uncorr)],
    by = c("aa_seq", "phenotype"))[!is.na(fitness)]
  pred_dt[, fitness_pred := mean]
  plot_dt <- pred_dt[!is.na(fitness_pred)&phenotype==phenotypen,.(observed_fitness = fitness, predicted_fitness = fitness_pred)]
  ggplot2::ggplot(plot_dt,ggplot2::aes(observed_fitness, predicted_fitness)) +
    ggplot2::stat_binhex(bins = 50, size = 0, color = "lightgrey") +
    ggplot2::scale_fill_gradientn(colours = c("white", "black"), trans = "log10") +
    # ggplot2::geom_point(alpha = 1/10) +
    ggplot2::xlab("Observed fitness (Rep. 3)") +
    ggplot2::ylab("Predicted fitness (Reps. 1&2)") +
    ggplot2::geom_text(data = plot_dt[,.(label = paste("r = ", round(cor(observed_fitness, predicted_fitness, use = "pairwise.complete"), 2), sep="")),],
                       x=-1,y=0.7,
                       ggplot2::aes(label=label, x=-Inf, y=Inf, hjust = 0, vjust = 1)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme_classic()+
    ggplot2::geom_abline(linetype="dashed")
}
rm(list = c("Plot_pre_12_ob_3_fitness_cor"))
Plot_pre_12_ob_3_fitness_cor(prediction="/Users/cweng/Desktop/Chenchun_KRAS_paper/revision/fullRAF/Mochi/predicted_phenotypes_all.txt",
                             observed="./DATA/predicted_phenotypes_all.txt",
                             phenotypen=1)





  
