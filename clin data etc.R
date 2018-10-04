


####mutational load by ER status
my_comparisons<-list(c("Negative", "Positive"))

g<-ggplot(comp[which(comp$ER.Status.By.IHC=="Positive" | comp$ER.Status.By.IHC=="Negative"),], aes(ER.Status.By.IHC, log(Mutation.Count)))
g+geom_point(position = "jitter",alpha=0.5)+
  stat_summary(fun.y = mean, geom = "bar", na.rm = TRUE, col="black",alpha=0.5)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.1, col="red")+
  geom_signif(comparisons = my_comparisons, 
              map_signif_level=TRUE, test="t.test", y_position = c(8), textsize = 8)+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  labs(x="ER status by IHC", y="Mutation count (log scale)")
table(comp$ER.Status.By.IHC [which(comp$ER.Status.By.IHC=="Positive" | comp$ER.Status.By.IHC=="Negative")])
###Mutational load relating to STAT1, IRF1, CD274 expression
which(comp$ER.Status.By.IHC=="Positive" | comp$ER.Status.By.IHC=="Negative")

gmutloadSTAT1<-ggplot(data = comp, aes(x=STAT1exp, y=log(Mutation.Count)))
my_comparisons<-list(c("Low", "Medium"), c("Low", "High"), c("Medium", "High"))
formula=y~x
par(mfrow=c(2,2))
gSTAT1<-gmutloadSTAT1+
  geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  geom_signif(comparisons = my_comparisons, 
              map_signif_level=TRUE, test="t.test", y_position = c(8, 10, 12), textsize = 8)+
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
ggplot(comp, aes(STAT1, log(Mutation.Count)))+geom_point()+
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x)+
  stat_poly_eq(formula=formula, aes(label = paste(..rr.label..)), parse=TRUE, size=8, label.x.npc = "right",
               label.y.npc = "bottom")+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  stat_fit_glance(method = 'lm',
                  method.args = list(formula = formula),
                  geom = 'text',
                  label.y.npc = "top",
                  label.x.npc = "right",
                  aes(label = paste("P = ", signif(..p.value.., digits = 2), sep = "")), size=8)

gmutloadIRF1<-ggplot(data = comp, aes(x=IRF1exp, y=log(Mutation.Count)))
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
  geom_signif(comparisons = my_comparisons, 
              map_signif_level=TRUE, test="t.test", y_position = c(8, 10, 12), textsize = 8)+
  labs(x='IRF1 expression', y='Mutational load (log scale)')+
  ylim(0,13)
gIRF1

gmutloadCD274<-ggplot(data = comp, aes(x=CD274exp, y=log(Mutation.Count)))
gCD274<-gmutloadCD274+
  geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  aes(x = fct_inorder(CD274exp))+
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14), 
        axis.text.y = element_text(size=14), 
        axis.title.x = element_text(size=18,
                                    margin = margin(t = 10, r = 0, b = 10, l = 0)), axis.title.y = element_text(size=18),
        legend.title = element_text(size=18),
        legend.text = element_text(size = 18))+
  labs(x='CD274 expression', y='Mutational load (log scale)')
grid.arrange(gSTAT1, gIRF1, gCD274, nrow = 1)



###Clinical
font.x=20, font.y=20, font.tickslab=20, font.legend= 15, font.pval=15
comp$ML<-fct_relevel(comp$ML, f)
library(survival)
library(survminer)
g <- survfit( Surv(Overall.Survival..Months., Status) ~ ML, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gML<-ggsurvplot(g, data=comp[which(comp$ER.Status.By.IHC=="Positive"),], pval = TRUE, pval.method = TRUE,
                legend.title="Tumour mutational load", legend.labs=c("Low", "Medium", "High"),
                palette = c("green", "blue", "red"), 
                xlab="Time (months)", ylab="Overall survival",
                font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=20)
gML
g <- survfit( Surv(Overall.Survival..Months., Status) ~ IRF1exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gIRF1<-ggsurvplot(g, data=comp[which(comp$ER.Status.By.IHC=="Positive"),], pval = TRUE, pval.method = TRUE,
                  legend.title="IRF1 expression", legend.labs=c("Low", "Medium", "High"),
                  palette = c("green", "blue", "red"), 
                  xlab="Time (months)", ylab="Overall survival",
                  font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=20)
gIRF1
?survdiff
g <- survfit( Surv(Overall.Survival..Months., Status) ~ STAT1exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gSTAT1<-ggsurvplot(g, data=comp[which(comp$ER.Status.By.IHC=="Positive"),], pval = TRUE, pval.method = TRUE,
                   legend.title="STAT1 expression", legend.labs=c("Low", "Medium", "High"),
           palette = c("green", "blue", "red"), 
           xlab="Time (months)", ylab="Overall survival",
           font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=20)
gSTAT1

g <- survfit( Surv(Overall.Survival..Months., Status) ~ CD274exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gCD274<-ggsurvplot(g, data=comp[which(comp$ER.Status.By.IHC=="Positive"),], pval = TRUE, pval.method = TRUE,
                   legend.title="CD274 expression", legend.labs=c("Low", "Medium", "High"),
                   palette = c("green", "blue", "red"), 
                   xlab="Time (months)", ylab="Overall survival",
                   font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=20)
gCD274
survdiff( Surv(Overall.Survival..Months., Status) ~ CD274exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
arrange_ggsurvplots(list(gML,gSTAT1,gCD274,gIRF1), nrow=2, ncol =  print=TRUE)


###DFS
g <- survfit( Surv(Disease.Free..Months., DStatus) ~ ML,
              data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
survdiff(Surv(Disease.Free..Months., DStatus) ~ ML,
         data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gML<-ggsurvplot(g, comp[which(comp$ER.Status.By.IHC=="Positive"),], pval=TRUE, pval.method = TRUE, 
          legend.title="Mutational load", legend.labs=c("Low", "Medium", "High"),
           palette = c("green", "blue", "red"),
           xlab="Time (months)", ylab="Progression free survival",
           font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=20)
gML
##DFS based on pd-L1
g <- survfit( Surv(Disease.Free..Months., DStatus) ~ CD274exp, 
              data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gCD274<-ggsurvplot(g, comp[which(comp$ER.Status.By.IHC=="Positive"),], pval=TRUE, pval.method = TRUE, 
                   legend.title="CD274 expression", legend.labs=c("Low", "Medium", "High"),
                   palette = c("green", "blue", "red"),
                   xlab="Time (months)", ylab="Progression free survival",
                   font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=20)
gCD274

##DFS based on IRF1
g <- survfit( Surv(Disease.Free..Months., DStatus) ~ IRF1exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gIRF1<-ggsurvplot(g, comp[which(comp$ER.Status.By.IHC=="Positive"),], pval=TRUE, pval.method = TRUE, pval.size=8, 
                  legend.title="IRF1 expression", legend.labs=c("Low", "Medium", "High"),
                  palette = c("green", "blue", "red"),
                  xlab="Time (months)", ylab="Progression free survival",
                  font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=15)
gIRF1
##DFS based on STAT1
g <- survfit( Surv(Disease.Free..Months., DStatus) ~ STAT1exp, data=comp[which(comp$ER.Status.By.IHC=="Positive"),])
gSTAT1<-ggsurvplot(g, comp[which(comp$ER.Status.By.IHC=="Positive"),], pval=TRUE, pval.method = TRUE, pval.size=8, 
                   legend.title="STAT1 expression", legend.labs=c("Low", "Medium", "High"), legend.position="left",
                   palette = c("green", "blue", "red"),
                   xlab="Time (months)", ylab="Progression free survival",
                   font.caption=20, font.x=20, font.y=20, font.tickslab=20, font.legend=15)
gSTAT1
arrange_ggsurvplots(list(gML,gSTAT1,gCD274,gIRF1), nrow=2, ncol=2, print=TRUE)
##########
summarySTAT<-comp %>% group_by(ML) %>% 
  mutate(mean=mean(`STAT1`), SE=std(`STAT1`), median = median(STAT1))
gSTAT1<-ggplot(data = comp, aes(x=ML, y=STAT1))
gSTAT1<-gSTAT1+geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  aes(x = fct_inorder(ML))+
  theme_bw()+
  labs(x='Mutational load', y='STAT1 expression')
gSTAT1
summaryIRF1<-comp %>% group_by(ML) %>% 
  mutate(mean=mean(`IRF1`), SE=std(`IRF1`), median = median(IRF1))
gIRF1<-ggplot(data = comp, aes(x=ML, y=IRF1))
gIRF1<-gIRF1+geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  aes(x = fct_inorder(ML))+
  theme_bw()+
  labs(x='Mutational load', y='IRF1 expression')
gIRF1
summaryCD274<-comp %>% group_by(ML) %>% 
  mutate(mean=mean(`CD274`), SE=std(`CD274`), median = median(CD274))
gCD274<-ggplot(data = comp, aes(x=ML, y=CD274))
gCD274<-gCD274+geom_point(stat='identity', position='jitter')+
  geom_boxplot(alpha=0.5)+
  aes(x = fct_inorder(ML))+
  theme_bw()+
  labs(x='Mutational load', y='CD274 expression')
gCD274


