#-----------------Creating plots from raw results of SPORTS1.0-----------------

library(tidyverse);library(ggpmisc);library(gridExtra);library(wesanderson)
mysum <- read_delim("13_S13_R1_001_summary.txt","\t")
mysum %>% 
  select(Class,Reads) %>% 
  filter(str_detect(Class,"Match_Genome|Clean")) %>% 
  separate(col=Class,into = c("Database","sRNA class"),sep = "-") ->myTb
myTb %>% 
  mutate(Database= str_replace(Database,"_Match_Genome","")) %>%
  mutate(Database= str_replace(Database,"Unannotated","Not annotated")) %>% 
           filter(!str_detect(Database,"Match_Genome|Clean"),
                  is.na(`sRNA class`)|!str_detect(`sRNA class`,"tRNA_5_end") ) %>%
  group_by(Database) %>% 
  summarise(reads=sum(Reads)) ->dat1

myTb$`sRNA class` <-sub("Match_Genome","",myTb$`sRNA class`)
myTb %>% 
  mutate( `sRNA class`= str_replace(`sRNA class`,"_Match_Genome","")) %>% 
  mutate(`sRNA class`=str_replace_all(`sRNA class` ,"_"," ")) %>% 
  group_by(`sRNA class`) %>% 
  summarise(reads=sum(Reads)) %>% drop_na()->dat2


Title_name <- c("testis")
p=ggplot(dat, aes(fill=`sRNA class`, y=reads/10*6, x=`sRNA class`)) + 
  geom_bar(stat="identity",position = "dodge")
p+scale_fill_manual(values = pal1)
  theme_minimal()+theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
  ggtitle(paste0("sRNA classes ",Title_name)) + 
  theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0.5))+
  scale_y_continuous("Raw reads",labels = scales::comma)+coord_flip()->fplot

write_rds(fplot,)

#---------------------create pallete--------------------------------------------
pal1 <- c("#9986A5", "#79402E", "#CCBA72", "#FD6467",
          "#5B1A18", "#D9D0D3", "#8D8680")
#-----------------Creating  txt summaries and plots-----------------------------
setwd("/home/rstudio/Documents/SPORTS_results/Result_furlan_with_DBs/")#change it#
filesnames <- list.files(pattern = "_summary.txt", recursive = TRUE)
dir.create("/home/rstudio/Documents/SPORTS_results/summary_results_furlan_with_DBs/")
#--------------name the folders-------------------------------------------------
dirs1 <- c("my_summaries","txt_dist", "sports_plots")
folders <- c("rRNAmap","rsRNA","tRNA","smallRNA")
#--------------create the folders-----------------------------------------------
fileout <- c("summary_results_furlan_with_DBs")# change it ##

lapply(paste0("/home/rstudio/Documents/SPORTS_results/",
              fileout,"/",dirs1),dir.create)

lapply(paste0("/home/rstudio/Documents/SPORTS_results/",
              fileout,"/sports_plots/",folders),dir.create)

for (filename in filesnames){
  filesnames1 <- gsub("_summary.txt","",filename)
  filesnames1 <- basename(filesnames1)
  mysum<- read_delim(filename,"\t")
  mysum %>% 
    select(Class,Reads) %>% 
    filter(str_detect(Class,"Match_Genome|Clean")) %>% 
    separate(col=Class,into = c("Database","sRNA class"),sep = "-") ->myTb
  write_delim(myTb,path = paste0("/home/rstudio/Documents/SPORTS_results/"
                                 ,fileout,
                                 "/my_summaries/",filesnames1,
                                 "_summary.txt"),delim = "\t")
  #--------------clean data for sRNA classes------------------------------------
  myTb %>% 
    mutate( `sRNA class`= str_replace(`sRNA class`,"_Match_Genome","")) %>% 
    mutate(`sRNA class`=str_replace_all(`sRNA class` ,"_"," ")) %>% 
    group_by(`sRNA class`) %>% 
    summarise(reads=sum(Reads)) %>% drop_na()->dat2
  #--------------ready the plot for sRNA classes--------------------------------
  p1=ggplot(dat2, aes(fill=`sRNA class`, 
                      y=reads/sum(dat2$reads), x=`sRNA class`)) + 
    geom_bar(stat="identity",position = "dodge")
  p1+
    theme_minimal()+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
    ggtitle(paste0("sRNA classes of sample:",filesnames1)) + 
    theme(plot.title = element_text(color="#666666",
                                    face="bold", size=18, 
                                    hjust=0.5),legend.position = "none")+
    geom_text(aes(label=scales::percent(reads/sum(dat2$reads)),
                  y=reads/sum(dat2$reads)),hjust=-.2)+
    scale_y_continuous("Relative frequencies of mapped reads",
                       labels = scales::percent,breaks = seq(0,1,by=0.05))+
    coord_flip()->fplot
  
  #-------save the code plot----------------------------------------------------
  write_rds(fplot,paste0("/home/rstudio/Documents/SPORTS_results/",fileout,
                         "/my_summaries/",filesnames1,"_plot_classes.rds"))
  #-------write the plot--------------------------------------------------------
  ggsave(filename = paste0(filesnames1,"_plot_classes.png"),plot=fplot,
         device = "png",
         path = paste0("/home/rstudio/Documents/SPORTS_results/",fileout,
         "/my_summaries/"),
         width= 180,
         height = 98,2,
         units = "mm",
         dpi = "retina")
  #---------------------clean data for Databases--------------------------------
  myTb %>% 
    mutate(Database= str_replace(Database,"_Match_Genome","")) %>%
    mutate(Database= str_replace(Database,"Unannotated","Not annotated")) %>% 
    filter(!str_detect(Database,"Match_Genome|Clean"),
           is.na(`sRNA class`)|!str_detect(`sRNA class`,"tRNA_5_end") ) %>%
    group_by(Database) %>% 
    summarise(reads=sum(Reads)) ->dat1
  #------------ready the plot for Databases-------------------------------------
  p2=ggplot(dat1, aes(fill=Database, y=reads/sum(dat1$reads), x=Database)) + 
    geom_bar(stat="identity",position = "dodge")
  p2+scale_fill_manual(values = pal1)+
    theme_minimal()+
    theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1),legend.position = "none")+
    ggtitle(paste0("Databases of ",filesnames1," sample")) + 
    theme(plot.title = element_text(color="#666666", face="bold", size=18, hjust=0.5))+
    geom_text(aes(label=scales::percent(reads/sum(dat1$reads)),
                  y=reads/sum(dat1$reads)),hjust=-.1)+
  scale_y_continuous("Relative frequencies of mapped reads",
                     labels = scales::percent,breaks = seq(0,1,by=0.05))+
    coord_flip()->oplot
  #---------------------save the code plot--------------------------------------
  write_rds(oplot,paste0("/home/rstudio/Documents/SPORTS_results/",fileout,
                         "/my_summaries/",filesnames1,"_plot_db.rds"))
  #---------------------write the plot------------------------------------------
  ggsave(filename = paste0(filesnames1,"_plot_db.png"),plot=oplot,device = "png",
         path = paste0("/home/rstudio/Documents/SPORTS_results/",fileout,
                       "/my_summaries/"),
         width= 180,
         height = 98,2,
         units = "mm",dpi = "retina")
} 



#-----copy all plot files in the other folder----------------------------------

filesnames_ld <- list.files(pattern = "_length_distribution.txt", recursive = TRUE)
filesnames_rRNAmap <- list.files(pattern = "_rRNA_mapping.pdf", recursive = TRUE)
filesnames_rsRNA<- list.files(pattern = "_rsRNA_distribution.pdf", recursive = TRUE)
filesnames_tRNA<- list.files(pattern = "_tRNA_mapping.pdf", recursive = TRUE)
filesnames_sncRNA<- list.files(pattern = "_sncRNA_distribution.pdf", recursive = TRUE)
input <- getwd()

#-------copy them--------------
#---------------------txt_dist
outfol <- paste0("/home/rstudio/Documents/SPORTS_results/",
                 fileout,"/txt_dist/")

file.copy(file.path(input,filesnames_ld),outfol, overwrite = "TRUE",
          recursive=TRUE)
#---------------------rRNAmap
outfol <- paste0("/home/rstudio/Documents/SPORTS_results/",
                 fileout,"/sports_plots/rRNAmap/")
file.copy(file.path(input,filesnames_rRNAmap),outfol)
#---------------------rsRNA
outfol <- paste0("/home/rstudio/Documents/SPORTS_results/",
                 fileout,"/sports_plots/rsRNA/")
file.copy(file.path(input,filesnames_rsRNA),outfol)

#---------------------tRNA
outfol <- paste0("/home/rstudio/Documents/SPORTS_results/",
                 fileout,"/sports_plots/tRNA/")
file.copy(file.path(input,filesnames_tRNA),outfol)

#---------------------smallRNA
outfol <- paste0("/home/rstudio/Documents/SPORTS_results/",
                 fileout,"/sports_plots/smallRNA/")
file.copy(file.path(input,filesnames_sncRNA),outfol)

#

#----------------------------Heatmap-------------------------------------------
library(ComplexHeatmap);library(tidyverse);library(edgeR);library(circlize)
# read first
setwd("/home/rstudio/Documents/new_analysis_workflow/results_R/results_with_DBs/Cell_lines_Testis_exp/")
DT<- read_delim(file = "all_mat_Cell_testis_noSK_norm_TMMwzp_logcpm_counts.txt",delim = "\t")
DT %>% column_to_rownames("smallRNA") -> DT

#DT %>% select(-Human_Testis) -> DTn
scaled_mat <- t(scale(t(as.matrix(DT))))
mat <- as.matrix(DT)
hist(scaled_mat)
df <- data.frame(cell_line = colnames(scaled_mat))
ha <- HeatmapAnnotation(df)
f <- colorRamp2(c(-1,median(scaled_mat),1),c("blue","white","yellow"))



Heatmap(scaled_mat,top_annotation = ha,
        col=f,show_row_dend = FALSE,clustering_distance_rows = "spearman",
        clustering_method_rows="ward.D2",clustering_distance_columns = "spearman",
        clustering_method_columns ="ward.D2",
        show_row_names = FALSE,name= "z-score equivalent", column_title = 
          "Heatmap of 4 piRNAS, only in cell lines",column_title_side="top")



# checking old new analysis files different versions of dbs--------
library(tidyverse);library(data.table)
setwd("~/Documents/piRNAs/complete_files_small_RNAs/")
fname <- list.files(pattern = "smallRNA_TOTAL")
old <- fread(fname[1],select = c(1,2,3,4,6)) %>% as_tibble()

fasta_old <- readDNAStringSet("/home/rstudio/Documents/SPORTS_Db/Homo_sapiens/piRBase/piR_human.fa")
fasta <- readDNAStringSet("/home/rstudio/Documents/SPORTS_Db/Homo_sapiens1/piRBase/piR_human.fa")

new <- DTRNA %>% 
  as_tibble() %>% 
  filter(str_detect(FileName,sub("smallRNA_","",sub(".txt","",fname[1]))))

old %>% filter(V6 %in% new$V6,str_detect(V4,"piR")) %>% arrange(desc(V6)) %>% view
new %>% filter(V6 %in% old$V6,str_detect(V4,"piR")) %>% arrange(desc(V6)) %>% view
