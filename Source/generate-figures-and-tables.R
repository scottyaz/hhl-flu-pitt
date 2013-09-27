## This file generates all the main findings for the manuscript
## Household Transmission of Influenza A and B in a School-based Study of Non-Pharmaceutical Interventions published in the Journal Epidemics

## Please contact aazman@jhsph.edu if you have any issues running the code
## we have posted this code knowing that it is not perfectly commented but believe that imperfectly commented code is better than no code.
## NOTE: this code is released under GPLv3 (http://www.gnu.org/licenses/gpl-3.0.html)

#########################
## Load some data files##
#########################

new.prim <- dget("Data/new_prim.rda")
prim.3 <- dget("Data/prim_3.rda") #individual level data for each primary case
new.prim <- dget("Data/new_prim.rda") # same as prim but suppoosedly with more reliable primary symptom data
hhls.in.prim.3 <- dget("Data/hhls_in_prim_3.rda") # indiviudal level data for household members and primary cases for all those in prim.3
sec.ILI.NoTemp <- dget("Data/sec_ILI_NoTemp.rda") # individual-level data on secondary cases linked to primary cases who meet the clincal definition of flu without the temerature requirement
sec.ILI <- dget("Data/sec_ILI.rda") # individual-level data on secondary cases linked to primary cases who meet the clincal definition of flu with the temerature requirement

########################
## Load some functions##
########################

source("Source/pipp-functions.R")
source("Source/likelihood-functions.R")

#####################################################
## Makes sets of ids for different case definitions##
#####################################################
#sets of ids for use in all analyses.  This allows us to limit to sets with different ILI definitions and those where the child is the index case
sets <- matrix(NA,nrow=7,ncol=2)
colnames(sets) <- c("n.hhls","n.people")
rownames(sets) <- c("all_pcr_pos","sec_notemp","sec_temp","sec_notemp_sickfirst","sec_temp_sickfirst","all_notemp_sickfirst","all_temp_sickfirst")
#Set 1: All households with PCR+ Children
ids.1 <- prim.3$ID
sets[1,1] <- length(ids.1)
sets[1,2] <- nrow(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.1,])
#NO Temp ILI (only ids of those with secondary cases in hhl)
ids.2 <- unique(sec.ILI.NoTemp$ID)
sets[2,1] <- length(ids.2)
sets[2,2] <- nrow(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.2,])
#Temp ILI (only ids of those with secondary cases in hhl)
ids.3 <- unique(sec.ILI$ID)
sets[3,1] <- length(ids.3)
sets[3,2] <- nrow(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.3,])
#Now eliminating households where student was not first case
#with NoTemp 28 households (only ids of those with secondary cases in hhl)
bef.aft.NoTemp <- before.after(Temp=FALSE)
ids.4 <- rownames(bef.aft.NoTemp[which(bef.aft.NoTemp[,4] == 1)  ,])
sets[4,1] <- length(ids.4)
sets[4,2] <- nrow(hhls.in.prim.3[which(hhls.in.prim.3$Sick.Childs.ID %in% ids.4),])
#with temp 24 hhls (only ids of those with secondary cases in hhl)
bef.aft.Temp <- before.after(Temp=TRUE)
ids.5 <- rownames(bef.aft.Temp[bef.aft.Temp[,4] == 1,])
sets[5,1] <- length(ids.5)
sets[5,2] <- nrow(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.5,])
#with NoTemp 71 households less than temp due to households where someone was sick (notemp) before child being counted in Temp only
# e.g. a nother child is "sick"" before the index case and since we are using a loose ILI def, they are counted as ILI and therefore this household is excluded
# in the Temp scenario this first case would not count as an ILI
ids.6 <- rownames(bef.aft.NoTemp[which(bef.aft.NoTemp[,4] == 1 | bef.aft.NoTemp[,4] == 3)  ,])
#1 = first case in hhl
#3 = no secondary cases
sets[6,1] <- length(ids.6)
sets[6,2] <- nrow(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.6,])
#with temp 74
ids.7 <- rownames(bef.aft.Temp[which(bef.aft.Temp[,4] == 1 | bef.aft.Temp[,4] == 3) ,])
sets[7,1] <- length(ids.7)
sets[7,2] <- nrow(hhls.in.prim.3[which(hhls.in.prim.3$Sick.Childs.ID %in% ids.7),])
sets

##################
## Make Table 1 ##
##################
table.1 <- make.table.1(ids=ids.1)  # the tab.type and tab.int elements make up the table
table.1$tab.type
table.1$tab.int

##################
## Make Table 2 ##
##################

## get SARs (these are from Table 2 and Table S2)
get.SARs(ids.1,by="INT")
get.SARs(ids.1,by="TYPE")

get.SARs(ids.7,by="TYPE") # table S2
get.SARs(ids.7,by="INT") # table S2

## get SITPs for final model used (other model fits can be generated using supplement code) assuming a family size of 4
SITPs <- get.SITPs(n.boots=500)

##################
## Make Table 3 ##
##################
## here we show the final GEE fit
library(geepack)
#restrict to secondary cases
prim.rows <- which(hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)")
secs.gee <- hhls.in.prim.3[-c(prim.rows),]

# or restrict to those who are first cases only
prim.rows.ind <- which(hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & hhls.in.prim.3$Sick.Childs.ID %in% ids.4)
secs.gee.firstinhhl <- hhls.in.prim.3[-c(prim.rows.ind),]
gee.df <- make.gee.df(secs.gee)
#gee.df = MakeGeeDf(secs.gee.firstinhhl)
#lose 14 observations and one full household
ids.exc <- which(is.na(gee.df$SchoolOrDaycare0n1y) | is.na(gee.df$Black0n1y))

gee.fit.4 <- geeglm(SICK.NoTemp ~ HouseholdSize + Cont0Int1 +ImmunThisYear0n1y + Black0n1y + OnlyResid0n1y + Flu0A1B + as.factor(Age.cat) , data=gee.df[-c(ids.exc),], id = as.factor(ID),family=binomial("logit"),corstr="exchangeable")

sum.gee.4 <- summary(gee.fit.4)
lower <- exp(sum.gee.4$geese$mean$estimate - qnorm(.975)*sum.gee.4$geese$mean$san.se)
upper <- exp(sum.gee.4$geese$mean$estimate + qnorm(.975)*sum.gee.4$geese$mean$san.se)
estimates <- exp(sum.gee.4$geese$mean$estimate)
tab.paper <- cbind(round(estimates,2),paste("(",round(lower,2),",",round(upper,2),")",sep=""),round(sum.gee.4$geese$mean$p,3))

# print out the table
data.frame(tab.paper)[-1,]

###################
## Make Figure 1 ##
###################
# pipp graphs
library(RColorBrewer)
load("Data/hhltables.rda")
#flutabs <- make.finaloutbreak.tables()

hhl.table.A.NoTemp<-hhl.table.A.NoTemp[c(-8),c(-6,-7,-8)]
hhl.table.B.NoTemp<-hhl.table.B.NoTemp[c(-8),c(-6,-7,-8)]
hhl.table.Inter.NoTemp<-hhl.table.Inter.NoTemp[c(-8),c(-6,-7,-8)]
hhl.table.Cont.NoTemp<-hhl.table.Cont.NoTemp[c(-8),c(-6,-7,-8)]

max(c(colSums(hhl.table.A.NoTemp), colSums(hhl.table.B.NoTemp), colSums(hhl.table.Inter.NoTemp), colSums(hhl.table.Cont.NoTemp)))

pals<-(brewer.pal(11,"Paired"))
blues<-(brewer.pal(6,"Blues"))
greens<-(brewer.pal(6,"Greens"))
reds<-(brewer.pal(6,"Reds"))

m <- .35
drop <- .91
cx <- 1.5

pltfnc<-function(hhl.table.A.NoTemp=hhl.table.A.NoTemp, hhl.table.B.NoTemp=hhl.table.B.NoTemp, pals=pals, greens=greens, blues=blues, annotate=T, LEG=c("A", "B")){

    ifelse(annotate, lab<-c("A", "B", "C"), lab<-c("D", "E", "F"))

    df<-rbind(c(colSums(hhl.table.A.NoTemp, na.rm=T)), c(colSums(hhl.table.B.NoTemp)))
    barplot(df, beside=T, border=grey(0), col=c(greens[4], blues[4]), names.arg=rep("", 5), ylim=c(0,20), yaxt='n')
    legend("topright", legend=LEG, col=c(greens[4], blues[4]), pch=c(15,15), bty='n', pt.cex=1.5)
    if(annotate) mtext("Households", side=2, line=2, las=0, cex=.85)
    axis(2, labels=annotate)
    text(par("usr")[1]+.25,par("usr")[4]*drop, lab[1], pos=4, cex=cx)

    mat.A<-hhl.table.A.NoTemp
    mat.A<-apply(mat.A, 2, function(x) x/sum(x))
    x<-barplot(mat.A, col=greens, names.arg=rep("", 5), yaxt='n')
    for(j in 1:5){
        y<-mat.A[,j]
        y[y==0]<-NA
        y[!is.na(y)]<-cumsum(y[!is.na(y)])
        y<-y[!is.na(y)]
        tt<-mat.A[,j]!=0
        text(x[j],y+.03, names(mat.A[,j])[tt][-1], cex=.95)
    }
    text(x, 0.05, 0, cex=.95)
    text(par("usr")[1]+.25,par("usr")[4]*drop, lab[2], pos=4, cex=cx)
                                        #legend("topright", legend=c(1:6,"SAR"), col=c(greens[1:6], pals[8]), pch=c(rep(15, 6), 17), bty='n', xpd=NA)
    sar.A<-hhl.table.A.NoTemp
    temp<-rep(NA, ncol(sar.A))
    for(j in 1:ncol(sar.A)) temp[j]<-sum(sar.A[,j]*0:4)/(sum(sar.A[,j])*(j+1))
    lines(x,temp, type='b', pch=24, col=grey(0), bg=pals[8], lwd=.5, cex=1.25)
                                        #axis(4, line=-3)
    if(annotate) mtext("Proportion", side=2, line=2, las=0, cex=.85, at=0)
    axis(2, labels=annotate)
    mat.B<-hhl.table.B.NoTemp
    mat.B<-apply(mat.B, 2, function(x) x/sum(x))
	barplot(mat.B, col=blues, yaxt='n')
    for(j in 1:5){
        y<-mat.B[,j]
        y[y==0]<-NA
        y[!is.na(y)]<-cumsum(y[!is.na(y)])
        y<-y[!is.na(y)]
        tt<-mat.B[,j]!=0
        text(x[j],y+.03, names(mat.B[,j])[tt][-1], cex=.95)
    }
    axis(2, labels=annotate)
    text(x, 0.05, 0, cex=.95)
    text(par("usr")[1]+.25,par("usr")[4]*drop, lab[3], pos=4, cex=cx)
                                        #legend("topright", legend=c(1:6,"SAR"), col=c(blues[1:6], pals[8]), pch=c(rep(15, 6), 17), bty='n', xpd=NA)
    sar.B<-hhl.table.B.NoTemp
    temp<-rep(NA, ncol(sar.B))
    for(j in 1:ncol(sar.B)) temp[j]<-sum(sar.B[,j]*0:4)/(sum(sar.B[,j])*(j+1))
    lines(x,temp, type='b', pch=24, col=grey(0), bg=pals[8], lwd=.5, cex=1.25)
                                        #axis(4, line=-3)

}

pals<-(brewer.pal(11,"Paired"))
browns<-rev((brewer.pal(11,"BrBG"))[1:6])
teals<-(brewer.pal(11,"BrBG"))[6:11]
blues<-(brewer.pal(6,"Blues"))
greens<-(brewer.pal(6,"Greens"))

purples<-rev((brewer.pal(11,"PRGn"))[1:6])
greens2<-(brewer.pal(11,"PRGn"))[6:11]
# reds<-rev((brewer.pal(11,"RdBu"))[1:6])
blues2<-(brewer.pal(11,"RdBu"))[6:11]

quartz("main", 5.5,6)
layout(matrix(1:6, byrow=F, ncol=2))
par(oma=c(3,3,1,1), mar=c(m,m,m+.1,m+.25), mgp=c(2,.5,0), cex.axis=.85, tck=-.02, las=1)
pltfnc(hhl.table.A.NoTemp, hhl.table.B.NoTemp, pals, blues2, greens2, annotate=T)
pltfnc(hhl.table.Inter.NoTemp, hhl.table.Cont.NoTemp, pals, purples, reds, annotate=F, LEG=c("Intervention","Control"))
mtext("Household Size", side=1, line=1.5, las=0, cex=.85, outer=T)
