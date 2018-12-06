
library(copula)
library(hydroGOF)
library(fitdistrplus)

# ESTIMATE MARGINAL distribution PARAMETERS BY MLE METHODE
fit.d=fitdist(m.spi,"weibull")
fit.d
plot(fit.d)
cdfcomp(list(fit.d),horizontals = TRUE,xlab="Duration data")

fit.s=fitdist(s.spi,"weibull")
fit.s
plot(fit.s)
cdfcomp(list(fit.s),horizontals = TRUE,xlab="Severity data")

# plot duration vs severity
plot(m.spi,s.spi,pch=1,xlab="Droght duration(month)",ylab="Droght severity")

#transform to uniform by fitted distribution
d.hat=pweibull (m.spi, shape=fit.d$estimate[1], scale=fit.d$estimate[2])
s.hat=pweibull (s.spi, shape=fit.s$estimate[1], scale=fit.s$estimate[2])
plot(d.hat,s.hat,xlab="duration-hat",ylab="severity-hat",title("Transformed data"))
uniform=cbind(d.hat,s.hat)

# transform  to pesduo observation by emperical distribution
m=cbind(m.spi,s.spi)
k=pobs(m)
plot(k,ylim=c(0,1),xlim=c(0,1))

# find intial guess for archemidan copula (gumbel,clayton,frank,joe)
gum.para  =initOpt("Gumbel",  interval=FALSE,  u=k)
clay.para =initOpt("Clayton", interval=FALSE,  u=k)
farnk.para=initOpt("Frank",   interval=FALSE,  u=k)
joe.para=  initOpt("Joe",     interval=FALSE,  u=k)

# create Dummy copula
gum.cop =     gumbelCopula (  param=gum.para,   dim=2)
clay.cop =    claytonCopula(  param=clay.para,  dim=2)
frank.cop =   frankCopula(    param=farnk.para, dim=2)
glambos.cop=  galambosCopula( param=2)
plackett.cop= plackettCopula( param=1)
joe.cop=      joeCopula(      param=joe.para,   dim=2)

# create marginal of copula
my.margins = c("weibull", "weibull")
my.params = list(list( shape=fit.d$estimate[1], scale=fit.d$estimate[2]),
                 list( shape=fit.s$estimate[1], scale=fit.s$estimate[2]))

# create mvdc class copula
myBvd.gum =      mvdc(copula=gum.cop,      margins=my.margins, paramMargins=my.params)
myBvd.caly =     mvdc(copula=clay.cop,     margins=my.margins, paramMargins=my.params)
myBvd.farnk =    mvdc(copula=frank.cop,    margins=my.margins, paramMargins=my.params)
myBvd.galambos=  mvdc(copula=glambos.cop,  margins=my.margins, paramMargins=my.params)
myBvd.plackett=  mvdc(copula=plackett.cop, margins=my.margins, paramMargins=my.params)
myBvd.joe =      mvdc(copula=joe.cop,      margins=my.margins, paramMargins=my.params)

##############################################################################
# Estimate Custom Bivariate Distn by MLE######################################
##############################################################################
start.vals1 = c(fit.d$estimate[1],fit.d$estimate[2],fit.s$estimate[1],fit.s$estimate[2],gum.para)
names(start.vals1)=c("m1", "m2", "m3", "m4", "m5")
myBvd.fitMvdc1 = fitMvdc(m, myBvd.gum, start.vals1)
my.params.fit1=list(list(shape=myBvd.fitMvdc1@estimate[1],scale=myBvd.fitMvdc1@estimate[2]),
                    list(shape=myBvd.fitMvdc1@estimate[3],scale=myBvd.fitMvdc1@estimate[4]))
d1=pweibull(m.spi,shape=myBvd.fitMvdc1@estimate[1],scale=myBvd.fitMvdc1@estimate[2])
s1=pweibull(s.spi,shape=myBvd.fitMvdc1@estimate[3],scale=myBvd.fitMvdc1@estimate[4])
u1=cbind(d1,s1)

start.vals2 = c(fit.d$estimate[1],fit.d$estimate[2],fit.s$estimate[1],fit.s$estimate[2],clay.para)
names(start.vals2)=c("m1", "m2", "m3", "m4", "m5")
myBvd.fitMvdc2 = fitMvdc(m, myBvd.caly, start.vals2)
my.params.fit2=list(list(shape=myBvd.fitMvdc2@estimate[1],scale=myBvd.fitMvdc2@estimate[2]),
                    list(shape=myBvd.fitMvdc2@estimate[3],scale=myBvd.fitMvdc2@estimate[4]))
d2=pweibull(m.spi,shape=myBvd.fitMvdc2@estimate[1],scale=myBvd.fitMvdc2@estimate[2])
s2=pweibull(s.spi,shape=myBvd.fitMvdc2@estimate[3],scale=myBvd.fitMvdc2@estimate[4])
u2=cbind(d2,s2)

start.vals3 = c(fit.d$estimate[1],fit.d$estimate[2],fit.s$estimate[1],fit.s$estimate[2],farnk.para)
names(start.vals3)=c("m1", "m2", "m3", "m4", "m5")
myBvd.fitMvdc3 = fitMvdc(m, myBvd.farnk, start.vals3)
my.params.fit3=list(list(shape=myBvd.fitMvdc3@estimate[1],scale=myBvd.fitMvdc3@estimate[2]),
                    list(shape=myBvd.fitMvdc3@estimate[3],scale=myBvd.fitMvdc3@estimate[4]))
d3=pweibull(m.spi,shape=myBvd.fitMvdc3@estimate[1],scale=myBvd.fitMvdc3@estimate[2])
s3=pweibull(s.spi,shape=myBvd.fitMvdc3@estimate[3],scale=myBvd.fitMvdc3@estimate[4])
u3=cbind(d3,s3)

start.vals4 = c(fit.d$estimate[1],fit.d$estimate[2],fit.s$estimate[1],fit.s$estimate[2],2)
names(start.vals4)=c("m1", "m2", "m3", "m4", "m5")
myBvd.fitMvdc4 = fitMvdc(m, myBvd.galambos, start.vals4)
my.params.fit4=list(list(shape=myBvd.fitMvdc4@estimate[1],scale=myBvd.fitMvdc4@estimate[2]),
                    list(shape=myBvd.fitMvdc4@estimate[3],scale=myBvd.fitMvdc4@estimate[4]))
d4=pweibull(m.spi,shape=myBvd.fitMvdc4@estimate[1],scale=myBvd.fitMvdc4@estimate[2])
s4=pweibull(s.spi,shape=myBvd.fitMvdc4@estimate[3],scale=myBvd.fitMvdc4@estimate[4])
u4=cbind(d4,s4)

start.vals5 = c(fit.d$estimate[1],fit.d$estimate[2],fit.s$estimate[1],fit.s$estimate[2],60)
names(start.vals5)=c("m1", "m2", "m3", "m4", "m5")
myBvd.fitMvdc5 = fitMvdc(m, myBvd.plackett, start.vals5)
my.params.fit5=list(list(shape=myBvd.fitMvdc5@estimate[1],scale=myBvd.fitMvdc5@estimate[2]),
                    list(shape=myBvd.fitMvdc5@estimate[3],scale=myBvd.fitMvdc5@estimate[4]))
d5=pweibull(m.spi,shape=myBvd.fitMvdc5@estimate[1],scale=myBvd.fitMvdc5@estimate[2])
s5=pweibull(s.spi,shape=myBvd.fitMvdc5@estimate[3],scale=myBvd.fitMvdc5@estimate[4])
u5=cbind(d5,s5)

start.vals6 = c(fit.d$estimate[1],fit.d$estimate[2],fit.s$estimate[1],fit.s$estimate[2],joe.para)
names(start.vals6)=c("m1", "m2", "m3", "m4", "m5")
myBvd.fitMvdc6 = fitMvdc(m, myBvd.joe, start.vals6)
my.params.fit6=list(list(shape=myBvd.fitMvdc6@estimate[1],scale=myBvd.fitMvdc6@estimate[2]),
                    list(shape=myBvd.fitMvdc6@estimate[3],scale=myBvd.fitMvdc6@estimate[4]))
d6=pweibull(m.spi,shape=myBvd.fitMvdc6@estimate[1],scale=myBvd.fitMvdc6@estimate[2])
s6=pweibull(s.spi,shape=myBvd.fitMvdc6@estimate[3],scale=myBvd.fitMvdc6@estimate[4])
u6=cbind(d6,s6)

myBvd.fitMvdc1
myBvd.fitMvdc2
myBvd.fitMvdc3
myBvd.fitMvdc4
myBvd.fitMvdc5
myBvd.fitMvdc6

gum.cop.fit =    gumbelCopula   (param=myBvd.fitMvdc1@estimate[5],  dim=2)
clay.cop.fit =   claytonCopula  (param=myBvd.fitMvdc2@estimate[5],  dim=2)
frank.cop.fit =  frankCopula    (param=myBvd.fitMvdc3@estimate[5],  dim=2)
glambos.cop.fit= galambosCopula (param=myBvd.fitMvdc4@estimate[5])
plackett.cop.fit=plackettCopula (param=myBvd.fitMvdc5@estimate[5])
joe.cop.fit=     joeCopula      (param=myBvd.fitMvdc6@estimate[5],  dim=2)

myBvd.fit1 = mvdc(copula=gum.cop.fit,       margins=my.margins, paramMargins=my.params.fit1)
myBvd.fit2 = mvdc(copula=clay.cop.fit,      margins=my.margins, paramMargins=my.params.fit2)
myBvd.fit3 = mvdc(copula=frank.cop.fit,     margins=my.margins, paramMargins=my.params.fit3)
myBvd.fit4 = mvdc(copula=glambos.cop.fit,   margins=my.margins, paramMargins=my.params.fit4)
myBvd.fit5 = mvdc(copula=plackett.cop.fit,  margins=my.margins, paramMargins=my.params.fit5)
myBvd.fit6 = mvdc(copula=joe.cop.fit,       margins=my.margins, paramMargins=my.params.fit6)
 
(gof.test1 = gofCopula(myBvd.fit1@copula, u1 ,estim.method="ml"))
(gof.test2 = gofCopula(myBvd.fit2@copula, u2 ,estim.method="ml"))
(gof.test3 = gofCopula(myBvd.fit3@copula, u3 ,estim.method="ml"))
(gof.test4 = gofCopula(myBvd.fit4@copula, u4 ,estim.method="ml"))
(gof.test5 = gofCopula(myBvd.fit5@copula, u5 ,estim.method="ml"))
(gof.test6 = gofCopula(myBvd.fit6@copula, u6 ,estim.method="ml")) 