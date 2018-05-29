library(ggplot2); library(gridExtra); library(lsmeans)

wes<-read.csv("data.csv")

###
#Figure 3 in McNickle and Evans 2018, AoB Plants. 
##

A<- ggplot(wes, aes(x = Treatment, y = lnRRfruit)) + geom_jitter(height = 0, width = 0.1, aes(alpha=0.35)) +
	theme(legend.position="none") + ylab("lnRR Fruit") + xlab("") +
	labs(title = "(A)") + geom_smooth(aes(x = as.numeric(Treatment), y = lnRRfruit))
B<- ggplot(wes, aes(x = Treatment, y = lnRRshoot)) +  geom_jitter(height = 0, width = 0.1, aes(alpha=0.35)) + 
	#geom_point(aes(y=pred_lnRRshoot), color="red", size=4) + 
	theme(legend.position="none") + ylab("lnRR Shoot") + 
	xlab("") + labs(title = "(B)") + geom_smooth(aes(x = as.numeric(Treatment), y = lnRRshoot)) +
	geom_smooth(aes(x = as.numeric(Treatment), y = pred_lnRRshoot, colour="red"))
C<- ggplot(wes, aes(x = Treatment, y = lnRRroot))  + geom_jitter(height = 0, width = 0.1, aes(alpha=0.35)) +
	theme(legend.position="none") + ylab("lnRR Root") + xlab("Damage treatment") +
	labs(title = "(C)") + geom_smooth(aes(x = as.numeric(Treatment), y = lnRRroot))
p15<-p5 + geom_jitter(height = 0, width = 0.1, aes(alpha=0.35))

dev.new()
grid.arrange(A, B, C, ncol=1)

##################################
###STATS reported in McNickle and Evans 2018, AoB Plants. 
##################################
library(lme4); library(car); library(pbkrtest);

pred <- lmer(pred_lnRRshoot~damage + (1|Block), data=wes)
obs <- lmer(lnRRshoot~damage + (1|Block), data=wes)
Anova(obs, test="F", type="III", ddf = "Kenward-Roger")
anova(obs,pred)

fruit.mod <- lmer(lnRRfruit~as.factor(damage) + (1|Block), data=wes)
Anova(fruit.mod, test="F", type="III", ddf = "Kenward-Roger")
lsmeans(fruit.mod, list(pairwise ~ as.factor(damage)), adjust = "tukey")

shoot.mod <- lmer(lnRRshoot~damage + (1|Block), data=wes)
Anova(shoot.mod, test="F", type="III", ddf = "Kenward-Roger")

root.mod <- lmer(lnRRroot~damage + (1|Block), data=wes)
Anova(root.mod, test="F", type="III", ddf = "Kenward-Roger")
