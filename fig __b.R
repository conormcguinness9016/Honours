library(org.Hs.eg.db)
library(cgdsr)
library(tidyverse)
library(corrplot)
library(reshape2)
library(ggpmisc)
kegg <- org.Hs.egPATH2EG

mapped <- mappedkeys(kegg)
kegg2 <- as.list(kegg[mapped])
kegg3<-as.character(unlist(kegg2))


symbol <- as.list(org.Hs.egSYMBOL) ## symbol holds listed gene symbols 
Sym <- na.omit(unlist(symbol[match(kegg3, names(symbol))]))
genes<-c(as.vector(Sym), "CD274", "IRF1")
mycgds = CGDS("http://www.cbioportal.org/public-portal/")
BRCAstudy<-getCancerStudies(mycgds)[32,1]
BRCA<-getCaseLists(mycgds, BRCAstudy)[1,1]
BRCAmeasures<-getGeneticProfiles(mycgds, BRCAstudy)
BRCARNAseq<-BRCAmeasures[4,1]
BRCAmeasures
i <- 1
while(i <= length(Sym)){
  if(i == 1){
    z_pca <- matrix()
  }
  k <- i+499
  temp <- getProfileData(mycgds, Sym[i:k], BRCARNAseq,BRCA)
  message(paste("done for genes -->",i," to ",k))
  if(dim(temp)[1] != 0){
    z_pca <- cbind(z_pca,  temp)
  }
  i <- k+1
}
z_pca1 <- t(z_pca)
rem_pca <- as.numeric(which(apply(z_pca1,1,function(x) sum(is.na(x))) == ncol(z_pca1)))
z_pca1 <- z_pca1[-rem_pca,]
z_pca1 <-t(z_pca1)
z_pca1<-z_pca1[as.vector(which(abs(rowSums(z_pca1))>0)),]
BRCAexp<-data.frame(z_pca1)
BRCAmut<-BRCAmeasures[8,1]
i <- 1
while(i <= length(Sym)){
  if(i == 1){
    z_pca <- matrix()
  }
  k <- i+499
  temp <- getProfileData(mycgds, Sym[i:k], BRCAmut,BRCA)
  message(paste("done for genes -->",i," to ",k))
  if(dim(temp)[1] != 0){
    z_pca <- cbind(z_pca,  temp)
  }
  i <- k+1
}
temp <- getProfileData(mycgds, c("miRNA221"), BRCAmiRNA,BRCA)
z_pcamut<-as.matrix(z_pca)
z_pcamut[is.na(z_pcamut)]<-"0"
z_pcamut[z_pcamut=="NaN"]<-"0"
z_pcamut[z_pcamut!="0"]<-"1"
z_pcamut<-data.frame(z_pcamut)
rowSums(z_pcamut=="1")

BRCAexp$mutcount<-NA
BRCAexp$mutcount<-rowSums(z_pcamut=="1")[match(rownames(BRCAexp), rownames(z_pcamut))]
BRCAexp$mutcount

cor<-cor(BRCAexp, method="spearman")
melted <- melt(cor, id.var="CD274")

rm(cor)
CD274<-melted[which(melted$Var1=="CD274"),]
CD274$highlight<-"dull"

CD274$highlight[which(CD274$Var2=="STAT1"|CD274$Var2=="STAT2"|CD274$Var2=="STAT3"|CD274$Var2=="STAT4"|CD274$Var2=="STAT5A"|CD274$Var2=="STAT5B"|CD274$Var2=="STAT6")|CD274$Var2=="IRF1"]<-"cool"
CD274$gene<-" "
CD274$gene[which(CD274$Var2=="STAT1")]<-"STAT1"
CD274$gene[which(CD274$Var2=="STAT2")]<-"STAT2"
CD274$gene[which(CD274$Var2=="STAT3")]<-"STAT3"
CD274$gene[which(CD274$Var2=="STAT4")]<-"STAT4"
CD274$gene[which(CD274$Var2=="STAT5A")]<-"STAT5A"
CD274$gene[which(CD274$Var2=="STAT5B")]<-"STAT5B"
CD274$gene[which(CD274$Var2=="STAT6")]<-"STAT6"
CD274$gene[which(CD274$Var2=="IRF1")]<-"IRF1"
CD274<-CD274[which(CD274$Var2!="CD274"),]
CD274$Var2<-factor(CD274$Var2, levels=CD274$Var2[order(CD274$value)])
CD274$Var2<-fct_relevel(CD274$Var2, CD274$Var2[order(-CD274$value)])
ggplot(CD274, aes(Var2, value, label=gene))+
  geom_point()+
  geom_point(data=CD274[which(CD274$Var2=="STAT1"),], size=2, col="red")+
  geom_point(data=CD274[which(CD274$Var2=="IRF1"),], size=2, col="light blue")+
  geom_text(aes(label=gene),hjust=1.1, vjust=0)+
  guides(col=FALSE)+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  scale_color_manual(values=c("red", "black"))+
  labs(x="Gene", y="Spearman's correlation value
     with CD274 expression (r)")

g<-ggplot(BRCAexp, aes(STAT1, CD274))
g+geom_point()+
  geom_smooth(method = "lm", color="black", formula = y ~ x)+
  stat_poly_eq(formula=formula, aes(label = paste(..rr.label..)), parse=TRUE, size=8,
               label.x.npc = "right", label.y.npc = "bottom")+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text',
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")), size=8)+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  labs(x="STAT1 expression", y="CD274 expression")


