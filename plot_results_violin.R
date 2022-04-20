
library(argparser)
library(ggplot2)
library(data.table)
library(ggpubr)


quant <- function(x) {
        quantile(x,c(0.05))
    }

quant1 <- function(x) {
        quantile(x,c(0))
    }


prep_data <- function(infile,facto) {
    table <- read.table(infile)
    table_pr <-data.frame(matrix(ncol = 3, nrow = 0))

    for (i in 1:51){
        for (r in 1:nrow(table)){
            new_row=c(table[r,i],i-1)
            table_pr <- rbind(table_pr, new_row) 
          }
    }
    table_pr[,3]=facto
    colnames(table_pr)<-c("Clustering_rate", "SNP_threshold", "fact")
    return(table_pr)
}


plot_cl_rates <- function(master_table,SNP_t,order,filename,x_l) {

    list_tab = list()
    
    for (r in 1:nrow(master_table)){
        list_tab[[r]] <- prep_data(paste0(master_table$stem[r],".all_cl_r_concat"),master_table$fact[r])

    }

    tablef <- rbindlist(list_tab)

    tablev_s <- subset(tablef, tablef$SNP_threshold <= SNP_t)


    tablev_s$SNP_threshold <- as.factor(tablev_s$SNP_threshold)
    tablev_s$fact <- as.factor(tablev_s$fact)

    tablev_s$fact <- factor(tablev_s$fact, levels=order)


#p <- ggplot(tablev_s, aes(x=fact, y=Clustering_rate)) + ylim(0,1)+
#  geom_violin(trim=FALSE,draw_quantiles = c(0.025, 0.5, 0.975))+
#  facet_wrap(~ SNP_threshold, nrow=1)
#p


    p <- ggplot(tablev_s, aes(x=fact, y=Clustering_rate)) + ylim(0,1)+
        geom_boxplot(lwd=0.25,  outlier.shape = 19, outlier.size = 0.3,outlier.stroke = 0.3)+
        ylab("Clustering rate") + xlab(x_l)+ facet_wrap(~ SNP_threshold, nrow=1) +
        ggtitle("SNP threshold") + theme_bw()+
#        theme(plot.title = element_text(hjust = 0.5))+ 
        theme(plot.title = element_text(size=14, hjust=0.5),
            axis.title.x = element_text(size=15),
            axis.title.y = element_text(size=15), panel.border = element_rect(color = "black",fill = NA,size = 1),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle = 60, hjust = 1,size=7))  
       
    tab_res<-aggregate(tablef$Clustering_rate~ tablef$fact+ tablef$SNP_threshold, tablef, quant)
    tab_res1<-aggregate(tablef$Clustering_rate~ tablef$fact+ tablef$SNP_threshold, tablef, quant1)


    tr=c()
    tr1=c()
    tr2=c()
    tr3=c()
    
    for (o in order){
    
    tr = c(tr,min(subset(tab_res,tab_res[,1] == o &  tab_res[,3] >= 0.95)[,2]))
    tr1 = c(tr1,min(subset(tab_res,tab_res[,1] == o &  tab_res[,3] == 1)[,2]))  
    tr2 = c(tr2,min(subset(tab_res1,tab_res1[,1] == o &  tab_res1[,3] >= 0.95)[,2]))  
    tr3 = c(tr3,min(subset(tab_res1,tab_res1[,1] == o &  tab_res1[,3] == 1)[,2]))  

    }
    
    
   m_tab= list()
   q_list= list()
    for (r in 1:nrow(master_table)){
        m_tab[[r]]<- read.csv(paste0(master_table$stem[r],".all_ldist_concat.csv"), header=T)
        m_tab[[r]][,3]= master_table$fact[r]
        q_list[[r]]<- quantile(m_tab[[r]][,2],c(0.025,0.975))
    
    }
    
    tablef1 = rbindlist(m_tab)

    colnames(tablef1)<-c("ID","tbl", "fact")

    tablef1$fact <- as.factor(tablef1$fact)
    
    tablef1$fact <- factor(tablef1$fact, levels=order)

#violin plots
    p1<- ggplot(tablef1, aes(x=fact, y=tbl)) +geom_violin(adjust = 7, trim=FALSE) + xlab(x_l) +
        ylab("TBL (SNPs)") + theme_bw()+
        theme(plot.title = element_text(hjust = 0.5))+ 
        theme(plot.title = element_text(size=14, hjust=0.5),
            axis.title.x = element_text(size=15), axis.text.x = element_text(size=rel(1.5)),axis.text.y = element_text(size=rel(1.5)),
            axis.title.y = element_text(size=15), panel.border = element_rect(color = "black",fill = NA,size = 1),panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))  
            
    
    
    
    figure <- ggarrange(p, p1,labels = c("a", "b"),ncol = 1, nrow = 2)

    figure

    
    ggsave(paste0(filename,".png"))
    ggsave(paste0(filename,".pdf"))
    
    res <- list("order of simulation settings"= order, "95% clustered in 95% simulations SNP threshold" = tr, "95% clustered in 100% simulations SNP threshold" = tr2, "100% clustered in 95% simulations SNP threshold" = tr1, "100% clustered in 100% simulations SNP threshold" = tr3,"95% TBLs" = q_list)

    sink(paste0(filename,"_results_summary.txt"))
    print(res)
    sink()
    
    
    return(res)
    
    
    
    
    
    
}

                                                                         
                                          
parser <- arg_parser("compare different simulation settings, generate plots and a summary res table")
parser <- add_argument(parser,"-S", type="integer", default=5,help="plot clustering rates for all SNP threshold =< to this value")
parser <- add_argument(parser, "-f",nargs=1,help="stem for output files")
parser <- add_argument(parser,"-l",nargs=1,help="x labels for plots. If spaces in label use quotes(eg. -l \"x label\")")
parser <- add_argument(parser,"-o",nargs= 1,help="comma delimited list of labels for each file in the desired order, eg. -o 1,2,3,4. If there are spaces within label use quotes")
parser <- add_argument(parser,"-i", nargs= 1,help="colon delimited list of stem name for each simulation (including path) same order as in -o. Use quotes if spaces in path")
                                                                                  
                                          
args <- parse_args(parser)                                

                                      

order<-unlist(strsplit(args$o, ","))


m_li <-unlist(strsplit(args$i, ":"))
master_list= list()                       
    
for (l in 1:length(rle(m_li)$values)){
	master_list[[l]]=as.data.frame(matrix(c(m_li[l],order[l]),ncol=2))
}


master_table =rbindlist(master_list,idcol=F,use.names=F)

                   
                                                                                 

rownames(master_table)<-NULL
colnames(master_table)<-c("stem","fact")

plot_cl_rates(master_table,args$S,order,args$f,args$l)


