#Gordon G. McNickle, Joel S Brown, Doug Lynch 2012
#Department of Biological Sciences
#University of Illinois at CHicago
#mcnickle@uic.edu
rm(list=ls(all=TRUE))
library(scatterplot3d)
library(rootSolve)
graphics.off()

###############################################
#Parameters
###############################################
## PLANT
############################################### 
## Rn    = Root resource, N concentration uM
## Rc    = Shoot resource, CO2 concentration ppm
## u[2]  = root strategy
## u[1]  = shoot strategy
## crn   = root cost in N units 
## crc	 = root cost in C units
## csn	 = shoot cost in N units
## csc 	 = shoot cost in C units
## alph	 = represents C:N ratio. Or plant need for C relative to N. 
## beta  = 1-alph	
################################################
## HERBIVORE
################################################
## u[3]  = level of attack
## ch    = cost associated with u[3]
## a	 = encounter rate of herbivore and plant tissue
################################################


plant.game<-function(u) {
		with (as.list(params),	{
	#Plant growth under shoot herbivory	
	#NB: u[1] is shoot, u[2] is root
		Hs <- Rc*(1-exp(-u[1]))		#harvest of C
		Ps <- Hs - csc*u[1] - crc*u[2] - alph*u[2]*d  #Net profit C
		dPs.dr <- (-crc - alph*d)		#derivative of Ps by r
		dPs.ds <- Rc*exp(-u[1])-csc	#derivative of Ps by s	

		Hr <- Rn*(1-exp(-u[2]))		#harvest of N
		Pr <- Hr - csn*u[1] - crn*u[2] - beta*u[2]*d	#Net profit N
		dPr.dr <- Rn*exp(-u[2])-crn - beta*d	#derivative of Pr by r
		dPr.ds <- (-csn)		#derivative of Pr by s

	#Herbivore investment		
		dPh.duh	<- a*u[2]*d - ch
		
		#derivatives of whole plant G-function where, G=(Ps^alph)*(Pr^beta)
		dG.dur <- E*alph*(Ps^(alph-1)) * (Pr^beta)*dPs.ds +
			  E*beta*(Ps^alph) * (Pr^(beta-1))*dPr.ds
		dG.dus <- E*alph*(Ps^(alph-1)) * (Pr^beta)*dPs.dr +
			  E*beta*(Ps^alph) * (Pr^(beta-1))*dPr.dr
		return(c(dG.dur = dG.dur, dG.dus = dG.dus))
		})} 

plant.game.alone<-function(u) {		#without Herbivore
		with (as.list(params),	{
		#NB: u[1] is shoot, u[2] is root
		Hs <- Rc*(1-exp(-u[1]))		#harvest of C
		Ps <- Hs - csc*u[1] - crc*u[2]  #Net profit C
		dPs.dr <- (-crc)		#derivative of Ps by r
		dPs.ds <- Rc*exp(-u[1])-csc	#derivative of Ps by s	

		Hr <- Rn*(1-exp(-u[2]))		#harvest of N
		Pr <- Hr - csn*u[1] - crn*u[2]	#Net profit N
		dPr.dr <- Rn*exp(-u[2])-crn	#derivative of Pr by r
		dPr.ds <- (-csn)		#derivative of Pr by s
				
		#derivatives of whole plant G-function where, G=(Ps^alph)*(Pr^beta)
		dG.dur <- E*alph*(Ps^(alph-1)) * (Pr^beta)*dPs.ds +
			  E*beta*(Ps^alph) * (Pr^(beta-1))*dPr.ds
		dG.dus <- E*alph*(Ps^(alph-1)) * (Pr^beta)*dPs.dr +
			  E*beta*(Ps^alph) * (Pr^(beta-1))*dPr.dr
		return(c(dG.dur = dG.dur, dG.dus = dG.dus))
		})} 

#Loop to iterate ESS sover over different resource amounts. 

#Create empty vectors to store output
CO2<-as.numeric() 	#empty C vector
N<-as.numeric() 	#empty N vector
d.vec<-as.numeric() 
Root<-as.numeric()	#empty root vector
Shoot<-as.numeric()	#empty shoot vector
herb<-as.numeric()	#empty herbivore vector
yield<-as.numeric()
Root.alone<-as.numeric()	#empty root vector
Shoot.alone<-as.numeric()	#empty shoot vector
yield.alone<-as.numeric()
damage<-as.numeric()
gameoff.root<-as.numeric()

#loop parameters
d.int<-.1		#Interval to increase N by in loop
d.max<-35		#Max N in loop
		#Min N in loop
Rc<-5000		#Set Rc to Rc.min		
Rn.int<-5		#Set Rn to Rn.min
Rn.vec<-as.numeric()

alpha<-0.95		#alpha value where, beta=(1-alpha)


#Loop
for(i in 1:50) {
	Rn<-Rn.int+2*(i-1); d<--0.01
	while (d<d.max) {
		d<-d+d.int
		params<- c(E=0.1, Rc=Rc, Rn=Rn, crc=2, csc=3, crn=0.01, csn=0.03, alph=alpha, beta=(1-alpha), a=0.8, ch=1)
		d.vec<-c(d.vec,d)
		ESS<-multiroot(f = plant.game, start = c(.1,.1), maxiter=10000, positive=TRUE)$root
		ESS.alone<-multiroot(f = plant.game.alone, start = c(.1,.1), maxiter=10000, positive=TRUE)$root
		r.temp<-ESS[2]
		s.temp<-ESS[1]
		r.temp.alone<-ESS.alone[2]
		s.temp.alone<-ESS.alone[1]
		Root<-c(Root,r.temp)
		Shoot<-c(Shoot,s.temp)
		herb<-c(herb,ESS[3])
		Root.alone<-c(Root.alone,r.temp.alone)
		Shoot.alone<-c(Shoot.alone,s.temp.alone)
		Rn.vec<-c(Rn.vec,Rn)
		#YIELD CALCULATIONS
		a<-0.7
		ch<-3
		crc<-1.2
		csc<-3
		crn<-0.1
		csn<-0.1
		E=.1
		Hs <- Rc*(1-exp(-s.temp))		#harvest of C
		Ps <- Hs - csc*s.temp - crc*r.temp - alpha*r.temp*d 	#Net profit C
		Hr <- Rn*(1-exp(-r.temp))		#harvest of N
		Pr <- Hr - csn*s.temp - crn*r.temp - (1-alpha)*r.temp*d	#Net profit N	
		Hh <- r.temp*d
		y.temp <- E*(Ps^alpha)*(Pr^(1-alpha))
		yield<-c(yield, y.temp)
		Hs.alone <- Rc*(1-exp(-s.temp.alone*exp(-a*0)))		#harvest of C
		Ps.alone <- Hs.alone - csc*s.temp.alone - crc*r.temp.alone  	#Net profit C
		Hr.alone <- Rn*(1-exp(-r.temp.alone))		#harvest of N
		Pr.alone <- Hr.alone - csn*s.temp.alone - crn*r.temp.alone	#Net profit N	
		y.temp.alone <- E*(Ps.alone^alpha)*(Pr.alone^(1-alpha))
		yield.alone<-c(yield.alone, y.temp.alone)
		d.temp<-if(Hh==0) {0} else {(Hh)/(r.temp+Hh)}
		damage<-c(damage, d.temp)
		gameoff.root<-c(gameoff.root, (1-d.temp)*r.temp.alone)
		}
}

lnRR.root<-log(Root/Root.alone)
lnRR.shoot<-log(Shoot/Shoot.alone)
lnRR.yield<-log(yield/yield.alone)

pred_lnRR.root<-log(gameoff.root/Root.alone)

out<-data.frame(d.vec, Root, Shoot, herb, yield, Root.alone, Shoot.alone, yield.alone, 
	damage, lnRR.root, lnRR.shoot, lnRR.yield, Rn.vec, gameoff.root, pred_lnRR.root)
write.csv(out, "C:/Users/gmcnickle/Desktop/out.root.herb.csv") 

#plot results
library(ggplot2); library(gridExtra)

a1<-ggplot(out, aes(damage, lnRR.yield, colour=Rn.vec)) + geom_point() +
	theme(legend.position="none") + ylab("lnRR Fruit") + xlab("") +
	labs(title = "(D)   Below ground herbivory") 
b<-ggplot(out, aes(damage, lnRR.shoot, colour=Rn.vec)) + geom_point()+ 
	theme(legend.position="none") + ylab("lnRR Shoot") + xlab("") +
	labs(title = "(E)") 
c<-ggplot(out, aes(damage, lnRR.root, colour=Rn.vec)) + geom_point() +
	geom_point(aes(y=pred_lnRR.root), colour="red") + ylim(-1.5,0)+ 
	theme(legend.position="none") + ylab("lnRR Root") + xlab("Damage") +
	labs(title = "(F)") 
d<-ggplot(out, aes(d.vec, damage)) + geom_point() 


grid.arrange(a1,b,c,ncol=1)


