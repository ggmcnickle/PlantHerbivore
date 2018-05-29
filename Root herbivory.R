#Gordon G. McNickle, Wesley Evans
#Department of Botany and Plant Pathology
#Purdue University
#gmcnickle@purdue.edu

library(scatterplot3d)
library(rootSolve)


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
		Ps <- Hs - csc*u[1] - crc*u[2] - alph*u[2]*(1-exp(-a*u[3]))  #Net profit C
		dPs.dr <- (-crc - alph*(1-exp(-a*u[3])))		#derivative of Ps by r
		dPs.ds <- Rc*exp(-u[1])-csc	#derivative of Ps by s	

		Hr <- Rn*(1-exp(-u[2]))		#harvest of N
		Pr <- Hr - csn*u[1] - crn*u[2] - beta*u[2]*(1-exp(-a*u[3]))	#Net profit N
		dPr.dr <- Rn*exp(-u[2])-crn - beta*(1-exp(-a*u[3]))	#derivative of Pr by r
		dPr.ds <- (-csn)		#derivative of Pr by s

	#Herbivore investment		
		dPh.duh	<- a*u[2]*exp(-a*u[3]) - ch
		
		#derivatives of whole plant G-function where, G=(Ps^alph)*(Pr^beta)
		dG.dur <- alph*(Ps^(alph-1)) * (Pr^beta)*dPs.ds +
			  beta*(Ps^alph) * (Pr^(beta-1))*dPr.ds
		dG.dus <- alph*(Ps^(alph-1)) * (Pr^beta)*dPs.dr +
			  beta*(Ps^alph) * (Pr^(beta-1))*dPr.dr
		return(c(dG.dur = dG.dur, dG.dus = dG.dus, dPh.duh = dPh.duh))
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
		dG.dur <- alph*(Ps^(alph-1)) * (Pr^beta)*dPs.ds +
			  beta*(Ps^alph) * (Pr^(beta-1))*dPr.ds
		dG.dus <- alph*(Ps^(alph-1)) * (Pr^beta)*dPs.dr +
			  beta*(Ps^alph) * (Pr^(beta-1))*dPr.dr
		return(c(dG.dur = dG.dur, dG.dus = dG.dus))
		})} 

#Loop to iterate ESS sover over different resource amounts. 

#Create empty vectors to store output
CO2<-as.numeric() 	#empty C vector
N<-as.numeric() 	#empty N vector
Root<-as.numeric()	#empty root vector
Shoot<-as.numeric()	#empty shoot vector
herb<-as.numeric()	#empty herbivore vector
yield<-as.numeric()
Root.alone<-as.numeric()	#empty root vector
Shoot.alone<-as.numeric()	#empty shoot vector
yield.alone<-as.numeric()
damage<-as.numeric()

#loop parameters
Rn.int<-1		#Interval to increase N by in loop
Rn.max<-100		#Max N in loop
Rn.min<-5		#Min N in loop
Rc<-2000		#Set Rc to Rc.min		
Rn<-Rn.min		#Set Rn to Rn.min

alpha<-0.95		#alpha value where, beta=(1-alpha)

#Loop to iterate ESS sover over different resource amounts. 

#Create empty vectors to store output
CO2<-as.numeric() 	#empty C vector
N<-as.numeric() 	#empty N vector
Root<-as.numeric()	#empty root vector
Shoot<-as.numeric()	#empty shoot vector
herb<-as.numeric()	#empty herbivore vector
yield<-as.numeric()
Root.alone<-as.numeric()	#empty root vector
Shoot.alone<-as.numeric()	#empty shoot vector
yield.alone<-as.numeric()
damage<-as.numeric()

#loop parameters
Rn.int<-1		#Interval to increase N by in loop
Rn.max<-100		#Max N in loop
Rn.min<-10		#Min N in loop
Rc<-5000		#Set Rc to Rc.min		
Rn<-Rn.min		#Set Rn to Rn.min

alpha<-0.95		#alpha value where, beta=(1-alpha)


#Loop
	while (Rn<Rn.max) {
		Rn<-Rn+Rn.int
		params<- c(Rc=Rc, Rn=Rn, crc=2, csc=3, crn=0.01, csn=.03, alph=alpha, beta=(1-alpha), a=0.8, ch=1)
		N<-c(N,Rn)
		ESS<-multiroot(f = plant.game, start = c(.1,.1,.1), maxiter=10000, positive=TRUE)$root
		ESS.alone<-multiroot(f = plant.game.alone, start = c(.1,.1), maxiter=10000, positive=TRUE)$root
		r.temp<-(if (ESS[3]>0) {ESS[2]} else {multiroot(f = plant.game.alone, start = c(.1,.1), maxiter=10000, positive = TRUE)$root[2]})
		s.temp<-(if (ESS[3]>0) {ESS[1]} else {multiroot(f = plant.game.alone, start = c(.1,.1), maxiter=10000, positive = TRUE)$root[1]})
		r.temp.alone<-ESS.alone[2]
		s.temp.alone<-ESS.alone[1]
		Root<-c(Root,r.temp)
		Shoot<-c(Shoot,s.temp)
		herb<-c(herb,ESS[3])
		Root.alone<-c(Root.alone,r.temp.alone)
		Shoot.alone<-c(Shoot.alone,s.temp.alone)
		#YIELD CALCULATIONS
		a<-0.7
		ch<-3
		crc<-1.2
		csc<-3
		crn<-0.1
		csn<-0.1
		Hs <- Rc*(1-exp(-s.temp))		#harvest of C
		Ps <- Hs - csc*s.temp - crc*r.temp - alpha*r.temp*(1-exp(-a*ESS[3])) 	#Net profit C
		Hr <- Rn*(1-exp(-r.temp))		#harvest of N
		Pr <- Hr - csn*s.temp - crn*r.temp - (1-alpha)*r.temp*(1-exp(-a*ESS[3]))	#Net profit N	
		Hh <- r.temp*(1-exp(-a*ESS[3]))
		y.temp <- (Ps^alpha)*(Pr^(1-alpha))
		yield<-c(yield, y.temp)
		Hs.alone <- Rc*(1-exp(-s.temp.alone))		#harvest of C
		Ps.alone <- Hs.alone - csc*s.temp.alone - crc*r.temp.alone  	#Net profit C
		Hr.alone <- Rn*(1-exp(-r.temp.alone))		#harvest of N
		Pr.alone <- Hr.alone - csn*s.temp.alone - crn*r.temp.alone	#Net profit N	
		y.temp.alone <- (Ps.alone^alpha)*(Pr.alone^(1-alpha))
		yield.alone<-c(yield.alone, y.temp.alone)
		damage<-c(damage, if(Hh==0) {0} else {(Hh)/(r.temp+Hh)})
		}


lnRR.root<-log(Root/Root.alone)
lnRR.shoot<-log(Shoot/Shoot.alone)
lnRR.yield<-log(yield/yield.alone)

treat = rep("Root herbivore", length(Root))

out<-data.frame(treat, N, Root, Shoot, herb, yield, Root.alone, Shoot.alone, yield.alone, damage, lnRR.root, lnRR.shoot, lnRR.yield)
write.csv(out, "C:/Users/gmcnickle/Desktop/out.root.herb.csv") 

out=rbind(out, out2)

#plot results
library(ggplot2); library(gridExtra)

a<-ggplot(out, aes(N, lnRR.yield, colour=treat)) + geom_line() +
	ylab("lnRR Fruit") + xlab("") +
	labs(title = "(A)") + theme(legend.position=c(.69,.22), legend.title=element_blank())
b<-ggplot(out, aes(N, lnRR.shoot, colour=treat)) + geom_line() + 
	ylab("lnRR Shoot") + xlab("")+
	labs(title = "(B)")  + theme(legend.position="none")
c<-ggplot(out, aes(N, lnRR.root, colour=treat)) + geom_line()+
	ylab("lnRR Root") + xlab("N available") +
	labs(title = "(C)")  + theme(legend.position="none")
d<-ggplot(out, aes(N, damage, colour=treat)) + geom_line() +
	ylab("proportion damage \nby herbivore") + xlab("N available")+
	labs(title = "(D)")  + theme(legend.position="none") + ylim(0,0.6)

grid.arrange(a,b,c,d,ncol=2)


