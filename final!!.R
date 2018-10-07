library("tidyverse")
library(ggsignif)
library(data.table)
library(ggpmisc)
library(survival)
library(survminer)

std <- function(x) sd(x)/sqrt(length(x))
count<-read.delim("Mutation Count.txt") ##read in mutation count data for patients
mutcount<-na.omit(data.frame(count)) ##get rid of NAs
mutcount<-mutcount[which(mutcount$Mutation.Count!=0),] ##get rid of 0 mutational counts (log scale later)
exp<-read.delim("mRNA_exp.txt") ##mRNA data for JAK STAT genes from cbiolinks
CD274<-read.delim("CD274.txt") ##mRNA data for CD274
mutcount$STAT1<-exp[match(mutcount$Sample.ID,exp$SAMPLE_ID), 5] #bind STAT1 data to mutational load
mutcount$IRF1<-exp[match(mutcount$Sample.ID,exp$SAMPLE_ID), 7] #bind IRF1 data to mutational load
mutcount$CD274<-exp[match(mutcount$Sample.ID,CD274$SAMPLE_ID), 3]#bind CD274 data to mutational load
comp<-data.frame(mutcount)
comp<-comp[complete.cases(comp),]##bind all the rows
DT<-data.table(comp) ##data table for stratification
DT[order(STAT1),STAT1exp:=rep(c("Low", "Medium", "High"),each=320)] ##stratify patients based on STAT1 expression and assign labels
DT[order(IRF1),IRF1exp:=rep(c("Low", "Medium", "High"),each=320)] ##stratify patients based on IRF1 expression and assign labels
clin<-read.delim("CLIN.tsv") ##read in clinical data
comp<-data.frame(comp, clin[match(comp$Sample.ID,clin$Sample.ID), ], DT[match(comp$Sample.ID, DT$Sample.ID), 7:8]) ##add stratifications
f<-c("Low", "Medium", "High") 
comp$STAT1exp<-fct_relevel(as.factor(comp$STAT1exp), f) ##define levels of stratification
comp$IRF1exp<-fct_relevel(as.factor(comp$IRF1exp), f)
##DStatus
comp$DStatus<-NA
comp$DStatus[which(comp$Disease.Free.Status=="DiseaseFree")]<-0 ##assign binary disease status to patients
comp$DStatus[which(comp$Disease.Free.Status=="Recurred/Progressed")]<-1
##DFS based on STAT1 (Figure 3.12A)
g <- survfit( Surv(Disease.Free..Months., DStatus) ~ STAT1exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gSTAT1<-ggsurvplot(g, comp[which(comp$ER.Status.By.IHC=="Positive"),], pval=TRUE, pval.method = TRUE, pval.size=8, 
                   legend.title="STAT1 expression", legend.labs=c("Low", "Medium", "High"), legend.position="left",
                   palette = c("green", "blue", "red"),
                   xlab="Time (months)", ylab="Progression free survival",
                   font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=15)
gSTAT1

##DFS based on IRF1 (Figure 3.12B)
g <- survfit( Surv(Disease.Free..Months., DStatus) ~ IRF1exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gIRF1<-ggsurvplot(g, comp[which(comp$ER.Status.By.IHC=="Positive"),], pval=TRUE, pval.method = TRUE, pval.size=8, 
                  legend.title="IRF1 expression", legend.labs=c("Low", "Medium", "High"),
                  palette = c("green", "blue", "red"),
                  xlab="Time (months)", ylab="Progression free survival",
                  font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=15)
gIRF1
my_comparisons<-list(c("Low", "Medium"), c("Low", "High"), c("Medium", "High"))
formula=y~x
gmutloadSTAT1<-ggplot(data = comp[which(comp$ER.Status.By.IHC=="Positive"),], aes(x=STAT1exp, y=log(Mutation.Count))) 
##mutational load based on STAT1 expression (Figure 3.12D)
gSTAT1<-gmutloadSTAT1+
  geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, y_position = c(9,10.5,12), textsize=8, test = "t.test")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  labs(x='STAT1 expression', y='Mutational load (log scale)')+ylim(0,13)
gSTAT1
##mutational load based on IRF1 expression (Figure 3.12E)
gmutloadIRF1<-ggplot(data = comp[which(comp$ER.Status.By.IHC=="Positive"),], aes(x=IRF1exp, y=log(Mutation.Count)))
gIRF1<-gmutloadIRF1+
  geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  geom_signif(comparisons = my_comparisons, map_signif_level = TRUE, y_position = c(9,10.5,12), textsize=8, test = "t.test")+
  labs(x='IRF1 expression', y='Mutational load (log scale)')+
  ylim(0,13)
gIRF1


