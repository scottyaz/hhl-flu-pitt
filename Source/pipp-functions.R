
## Used to figure out temporal ordering of cases
## Helps consruct sets of of households where child is true index case
##' @param my.sec.ILI
##' @param my.sec.ILI.NoTemp
##' @param my.new.prim
##' @param my.prim.3
##' @param Temp
##' @return
##' @author Andrew Azman
before.after <- function(my.sec.ILI=sec.ILI,
                        my.sec.ILI.NoTemp=sec.ILI.NoTemp,
                        my.new.prim=new.prim,
                        my.prim.3=prim.3,
                        Temp=FALSE){
  if (Temp == TRUE){
    sec <- my.sec.ILI
  } else {
    sec <- my.sec.ILI.NoTemp
  }
  ids <- prim.3$ID
  z <- matrix(NA,nrow=length(ids),ncol=4)
  i <- 1
  for (j in ids){
    temp <- sec[sec$ID == j,]
    print(paste(j,nrow(temp)))
    z[i,3] <- my.new.prim$STSICKDT.y[my.prim.3$ID == j]
                                        #print(head(temp))
    if (nrow(temp) != 0){
    if (all(is.na(temp$STSICKDT))) {
      st <- NA
    } else {
      st <- min(temp$STSICKDT,na.rm=T)
    }
    z[i,1] <- st
    if (all(is.na(temp$SPSICKDT))) {
      sp <- NA
    } else {
      sp <- min(temp$SPSICKDT,na.rm=T)
    }
    z[i,2] <- sp
    if (is.na(st) & is.na(sp)) {
      z[i,4] <- 0  #0 for either missing information or someone else came first
    } else {
      z[i,4] <- ifelse(min(st,sp,na.rm=T) < my.new.prim$STSICKDT.y[my.prim.3$ID == j],0,1)
    }
  } else {
      z[i,4] <- 3 #3 means no secondary cases reported
    }
    i <- i+1
}
  rownames(z) <- ids
  return(z)
}

##Calculation of Secondary Attack Rates and exact binomial CIs
## ids - set of houshold ids to use (see sets object created above)
## ex.vac - exclude households where we are missing vaccination status?
## by = by INT means by intervention and by TYPE means by flutype
get.SARs <- function(ids=ids.1,ex.vac=FALSE,by="INT"){
  library("Hmisc")
  temp0 <- hhls.in.prim.3[which(hhls.in.prim.3$Sick.Childs.ID %in% ids),]
  temp0 <- temp0[-which(temp0$Relationship.to.Sick.Child == "Sick Child (Self)") ,]
  if(ex.vac){
    #going to exlude full households where we are missing any info on vaccination status
    exc.ids <- unique(temp0$Sick.Childs.ID[(which(is.na(temp0$Immunization.of.Family.Member)))])
    temp <- ddply(temp0[!temp0$Sick.Childs.ID %in% exc.ids ,],.(Sick.Childs.ID),summarise, "size"=Number.of.people.in.family[1]-1,"sick_notemp"=sum(SICK.NoTemp),"sick_temp"=sum(SICK.Temp),"int"=Intervention..y.n.[1],"type"=FluType[1])
  } else {
    temp <- ddply(temp0,.(Sick.Childs.ID),summarise, "size"=Number.of.people.in.family[1]-1,"sick_notemp"=sum(SICK.NoTemp),"sick_temp"=sum(SICK.Temp),"int"=Intervention..y.n.[1],"type"=FluType[1])
  }
  if (by=="INT"){
    sars <- ddply(temp,.(int),summarise,"sar_notemp"=sum(sick_notemp)/sum(size),"sar_temp"=sum(sick_temp)/sum(size),"num"=sum(size),"sick_notemp"=sum(sick_notemp),"sick_temp"=sum(sick_temp))
  } else if (by == "TYPE"){
    sars <- ddply(temp,.(type),summarise,"sar_notemp"=sum(sick_notemp)/sum(size),"sar_temp"=sum(sick_temp)/sum(size),"num"=sum(size),"sick_notemp"=sum(sick_notemp),"sick_temp"=sum(sick_temp))
} else {
    ##cut by int and type
    sars <- ddply(temp,.(type,int),summarise,"sar_notemp"=sum(sick_notemp)/sum(size),"sar_temp"=sum(sick_temp)/sum(size),"num"=sum(size),"sick_notemp"=sum(sick_notemp),"sick_temp"=sum(sick_temp))
    #stop("not using a valid 'by' type")
  }
  ci_temp <- binconf(sars$sick_temp,sars$num,method="exact") #may to use asymptotic CIs but doesn't really make a difference
  ci_notemp <- binconf(sars$sick_notemp,sars$num,method="exact")
  return(list("sars"=sars,"ci_temp"=ci_temp,"ci_notemp"=ci_notemp))
}


##' makes a list of final outbreak sizes for each comparison
##' need the id set from above in memory to run
##' @param first.cases.only
##' @return list of contingency tables
##' @author Andrew Azman
make.finaloutbreak.tables <- function(first.cases.only=FALSE){
   if (first.cases.only){
    table.base <- ddply(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.7,],.(Sick.Childs.ID),summarise, "hhlsize" = Number.of.people.in.family[1], "sec.sick.notemp" = (sum(SICK.NoTemp)-1), "sec.sick.temp" = (sum(SICK.Temp)-1), "flutype" = FluType[1], "int" = Intervention..y.n.[1])
    ##since there are a few households excluded when we use the NoTemp (due to missing info) we re-write over sick.notemps
    table.base2 <- ddply(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.6,],.(Sick.Childs.ID),summarise, "sec.sick.notemp" = (sum(SICK.NoTemp)-1), "sec.sick.temp" = (sum(SICK.Temp)-1), "flutype" = FluType[1], "int" = Intervention..y.n.[1])
    ##combine the dfs
    table.base.comb <- merge(table.base[,-3],table.base2[,c("Sick.Childs.ID","sec.sick.notemp")],by=c("Sick.Childs.ID"),all.x=T)
    table.base.comb$sec.sick.notemp <- factor(table.base.comb$sec.sick.notemp,levels=c(0,1,2,3,4,5))
    table.base.comb$sec.sick.temp <- factor(table.base.comb$sec.sick.temp,levels=c(0,1,2,3,4,5))
    table.base.comb$hhlsize <- factor(table.base.comb$hhlsize,levels=c(2,3,4,5,6))

    ##make tables
    flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,flutype))
    int.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int))
    flutype.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,flutype))
    int.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,int))
    int.flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int,flutype))
    flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,flutype))
    int.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int))
    flutype.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,flutype))
    int.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,int))
    int.flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int,flutype))


    flu.tables <- list(flutype.tables.NoTemp[,,1], #a
                       flutype.tables.NoTemp[,,2], #b
                       int.tables.NoTemp[,,1], #control
                       int.tables.NoTemp[,,2], #intervention
                       flutype.tables.Temp[,,1],
                       flutype.tables.Temp[,,2],
                       int.tables.Temp[,,1],
                       int.tables.Temp[,,2],
                       int.flutype.tables.NoTemp[,,1,1], #A cont
                       int.flutype.tables.NoTemp[,,2,1], #A int
                       int.flutype.tables.NoTemp[,,1,2], #B cont
                       int.flutype.tables.NoTemp[,,2,2]) #B int

  } else {
    table.base.comb <- ddply(hhls.in.prim.3,.(Sick.Childs.ID),summarise, "hhlsize" = Number.of.people.in.family[1], "sec.sick.notemp" = (sum(SICK.NoTemp)-1), "sec.sick.temp" = (sum(SICK.Temp)-1), "flutype" = FluType[1], "int" = Intervention..y.n.[1])
    ##combine the dfs
    table.base.comb$sec.sick.notemp <- factor(table.base.comb$sec.sick.notemp,levels=c(0:5))
    table.base.comb$sec.sick.temp <- factor(table.base.comb$sec.sick.temp,levels=c(0,1,2,3,4,5))
    table.base.comb$hhlsize <- factor(table.base.comb$hhlsize,levels=c(2,3,4,5,6))

    ##make tables
    flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,flutype))
    int.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int))
    flutype.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,flutype))
    int.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,int))
    int.flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int,flutype))
    flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,flutype))
    int.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int))
    flutype.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,flutype))
    int.tables.Temp <- with(table.base.comb,table(sec.sick.temp,hhlsize,int))
    int.flutype.tables.NoTemp <- with(table.base.comb,table(sec.sick.notemp,hhlsize,int,flutype))


    flu.tables <- list(flutype.tables.NoTemp[,,1], #a
                       flutype.tables.NoTemp[,,2], #b
                       int.tables.NoTemp[,,1], #control
                       int.tables.NoTemp[,,2], #intervention
                       flutype.tables.Temp[,,1],
                       flutype.tables.Temp[,,2], #
                       int.tables.Temp[,,1],
                       int.tables.Temp[,,2], #hhl size 5
                       int.flutype.tables.NoTemp[,,1,1], #A cont
                       int.flutype.tables.NoTemp[,,2,1], #A int
                       int.flutype.tables.NoTemp[,,1,2], #B cont
                       int.flutype.tables.NoTemp[,,2,2]) #B int
  }
   return(flu.tables)
}

##' This function estimates the SITPs for all different case definitions and cuts of the data
##' @param first.cases.only - TRUE if we want to only consider households where the school child was the first to report symptoms in the household
##' @return matrix with estimates and ci's
##' @author Andrew Azman
get.SITPs <- function(first.cases.only=FALSE,n.boots=500){
    load("../pipp/Rcode/hhltables.rda")
    mles <- matrix(NA,nrow=6,ncol=6)
    colnames(mles) <- c("SITP1","SITP2","ci1_U","ci1_L","ci2_U","ci2_L")
    rownames(mles) <- c("A/B.NoTemp","Cont/Int.NoTemp","A/B","Cont/Int","A.Cont/A.Int","B.Cont/B.Int")

    ##define the list of different tables that will be used:

    if (first.cases.only){
        flu.tables  <-  make.finaloutbreak.tables(first.cases.only=TRUE)
    } else {
        flu.tables  <-  make.finaloutbreak.tables(first.cases.only=FALSE)
    }

    firsts <- seq(1,11,by=2)
    seconds <- seq(2,12,by=2)

    for (i in 1:6){
        hhltab1 <- flu.tables[[firsts[i]]]
        hhltab2 <- flu.tables[[seconds[i]]]

        qs <- optim(c("q1"=.90,"q2"=.90,"alpha"=0),
                fn=lik3.2,
                table1=hhltab1,
                table2=hhltab2)$par

        mles[i,1:2] <- c(1-plogis(qs[1]/4^qs[3]),1-plogis(qs[2]/4^qs[3]))
        boots <- getBoots(hhltab1,hhltab2,n.boots)
        mles[i,3:4] <- apply(boots,2,function(x) quantile(x,c(0.025,0.975)))[,4]
        mles[i,5:6] <- apply(boots,2,function(x) quantile(x,c(0.025,0.975)))[,5]
    }
    return(mles)
}


## Functions to bootstrap contingency tables
##' @param table
##' @return
##' @author Andrew Azman
bootTable <- function(table){
    require(reshape)

    ## first take table and decompose into individual observations
    tmp <- melt(table)
    tmp <- tmp[which(!is.na(tmp[,3]) & tmp[,3] > 0 & !is.na(tmp[,1]) & !is.na(tmp[,2])),]
    colnames(tmp) <- c("sec.cases","hhl.size","count")
    tmp[,1] <- as.numeric(as.character(tmp[,1]))
    tmp[,2] <- as.numeric(as.character(tmp[,2]))
    tmp[,3] <- as.numeric(as.character(tmp[,3]))
    tmp <- as.matrix(tmp)
    ## first column is secondary cases, second column is household size, thrid is count
    ## now for each count we want to replicate an entry that many times
    ind.data <- c()
    for (i in 1:nrow(tmp)){
        ind.data <- rbind(ind.data,matrix(rep(tmp[i,1:2],tmp[i,3]),ncol=2,byrow=T))       }
    ## now resample
    new.ind.data <- ind.data[sample(nrow(ind.data),nrow(ind.data),replace=T),]

    ## now put it back into a table
    ## this is ineffecient but will work
    colnames <- as.numeric(colnames(table))
    if (any(is.na(colnames))) colnames <- colnames[-which(is.na(colnames))]
    rownames <- as.numeric(rownames(table))
    if (any(is.na(rownames))) rownames <- rownames[-which(is.na(rownames))]

    tmp.grid <- expand.grid(rownames,colnames)
    new.data <- c()
    for (j in 1:nrow(tmp.grid)){
        new.data[j] <- length(which(new.ind.data[,1] == tmp.grid[j,1] & new.ind.data[,2] == tmp.grid[j,2]))
    }

    out <- matrix(new.data,nrow=length(rownames))
    rownames(out) <- rownames
    colnames(out) <- colnames
    out
}

##' for each table we will resample final attack sizes then recompute the SITP
##' @param table1
##' @param table2
##' @return
singleBoot <- function(table1,table2){

    samp.tab1 <- bootTable(table1)
    samp.tab2 <- bootTable(table2)

    tmp <- optim(c("q1"=.90,"q2"=.90,"alpha"=0),fn=lik3.2,table1=samp.tab1,table2=samp.tab2)
    if (tmp$convergence != 1){
        c(tmp$par, SITP1= 1-plogis(tmp$par[1]/4^tmp$par[3]),SITP2=1-plogis(tmp$par[2]/4^tmp$par[3]))
    } else {
        rep(0,5)
    }
}

##' Runs n.boots bootstraps
##' @param tab1 contingency table 1
##' @param tab2 contengicy table 2
##' @param n.boots number of bootstraps
##' @return matrix of bootstrap
##' @author Andrew Azman
getBoots <- function(tab1,tab2,n.boots){
    rc <- array(dim=c(n.boots,5))
    for (i in 1:n.boots){
        if (i %% 50 == 0) cat(".")
    rc[i,] <- singleBoot(tab1,tab2)
    }
    rc
}

##' Makes table 1 for a given set of ids
##' Note: this relies on lexical scoping and should really only be called from the script so that everything that needs to be in memory is loaded correctly
##' @param ids
##' @return
##' @author Andrew Azman
make.table.1 <- function(ids){
    tab.int <- ddply(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids,],.(Intervention..y.n.),summarise,
                     "NoFluShot"=mean(Immunization.of.Family.Member == "No flu shot",na.rm=T),
                     "FluShotThisYear"=mean(Immunization.of.Family.Member == "This school year",na.rm=T),
                     "FluShotLastYear"=mean(Immunization.of.Family.Member == "Last school year",na.rm=T),
                     "PercentU15"=mean(Age.of.Family.Member <= 15,na.rm=T),
                     "PercentO15"=mean(Age.of.Family.Member > 15,na.rm=T),
                     "WorkFT"=mean(Activity.Out.of.Home.Fam.Member=="Work FT",na.rm=T),
                     "InSchoolorDayCare"=mean(Activity.Out.of.Home.Fam.Member=="School" | Activity.Out.of.Home.Fam.Member=="Daycare",na.rm=T),
                     "NotSchoolDayCareWorkFT"=mean(Activity.Out.of.Home.Fam.Member !="School" & Activity.Out.of.Home.Fam.Member !="Daycare",na.rm=T),
                     "Black"=mean(Race.of.Family.Member == "Black",na.rm=T),
                   "NotBlack"=mean(Race.of.Family.Member != "Black",na.rm=T),
                     "Chronic"=mean(Chronic.Illness..Fam.Member == "Yes",na.rm=T),
                     "NotChronic"=mean(Chronic.Illness..Fam.Member != "Yes",na.rm=T),
                     "FluA"=sum(FluType == "A" & Relationship.to.Sick.Child == "Sick Child (Self)",na.rm=T)/sum(Relationship.to.Sick.Child == "Sick Child (Self)",na.rm=T),
                   "MainChildAge"=mean(Age.of.Family.Member[Relationship.to.Sick.Child == "Sick Child (Self)"],na.rm=T))

  ##how many observations for each?
  tab.int.sizes <- ddply(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids,],.(Intervention..y.n.),summarise,
                         "NoFluShot"=sum(Immunization.of.Family.Member == "No flu shot",na.rm=T),
                         "FluShotThisYear"=sum(Immunization.of.Family.Member == "This school year",na.rm=T),
                         "FluShotLastYear"=sum(Immunization.of.Family.Member == "Last school year",na.rm=T),
                         "PercentU15"=sum(Age.of.Family.Member <= 15,na.rm=T),
                         "PercentO15"=sum(Age.of.Family.Member > 15,na.rm=T),
                         "WorkFT"=sum(Activity.Out.of.Home.Fam.Member=="Work FT",na.rm=T),
                         "InSchoolorDayCare"=sum(Activity.Out.of.Home.Fam.Member=="School" | Activity.Out.of.Home.Fam.Member=="Daycare",na.rm=T),
                         "NotSchoolDayCareWorkFT"=sum(Activity.Out.of.Home.Fam.Member !="School" & Activity.Out.of.Home.Fam.Member !="Daycare",na.rm=T),
                         "Black"=sum(Race.of.Family.Member == "Black",na.rm=T),
                         "NotBlack"=sum(Race.of.Family.Member != "Black",na.rm=T),
                         "Chronic"=sum(Chronic.Illness..Fam.Member == "Yes",na.rm=T),
                         "NotChronic"=sum(Chronic.Illness..Fam.Member != "Yes",na.rm=T),
                         "FluA"=sum(FluType == "A" & Relationship.to.Sick.Child == "Sick Child (Self)",na.rm=T ),
                         "MainChildAge"=sum(Relationship.to.Sick.Child == "Sick Child (Self)" & !is.na(Age.of.Family.Member) ,na.rm=T))

  num.cont <- ddply(hhls.in.prim.3[which(hhls.in.prim.3$Intervention..y.n. == "No" & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,],.(Sick.Childs.ID),function(x) x$Number.of.people.in.family[1])[,2]
  num.int <- ddply(hhls.in.prim.3[which(hhls.in.prim.3$Intervention..y.n. == "Yes" & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,],.(Sick.Childs.ID),function(x) x$Number.of.people.in.family[1])[,2]
tab.int <- cbind(tab.int,"Mean Household Size" = c(mean(num.cont),mean(num.int)))

                                        #prop tests
  i.tab <- vector("list",15)
  i.tab[[1]] <- test.NOFLUSHOT <- prop.test(tab.int.sizes[,2],c(sum(tab.int.sizes[1,2:4]),sum(tab.int.sizes[2,2:4])))
  i.tab[[2]] <- test.FluShotThisYear <-prop.test(tab.int.sizes[,3],c(sum(tab.int.sizes[1,2:4]),sum(tab.int.sizes[2,2:4])))
  i.tab[[3]] <- test.FluShotLastYear <-prop.test(tab.int.sizes[,4],c(sum(tab.int.sizes[1,2:4]),sum(tab.int.sizes[2,2:4])))
  i.tab[[4]] <- test.U15 <- prop.test(tab.int.sizes[,5],c(sum(tab.int.sizes[1,5:6]),sum(tab.int.sizes[2,5:6])))
  i.tab[[5]] <- test.O15 <- prop.test(tab.int.sizes[,6],c(sum(tab.int.sizes[1,5:6]),sum(tab.int.sizes[2,5:6])))
  i.tab[[6]] <- test.WorkFT <- prop.test(tab.int.sizes[,7],c(sum(tab.int.sizes[1,7:9]),sum(tab.int.sizes[2,7:9])))
  i.tab[[7]] <- test.INSCHOOLorDAYCARE <- prop.test(tab.int.sizes[,8],c(sum(tab.int.sizes[1,7:9]),sum(tab.int.sizes[2,7:9])))
  i.tab[[8]] <- test.NotSchoolDayCareWorkFT <- prop.test(tab.int.sizes[,9],c(sum(tab.int.sizes[1,7:9]),sum(tab.int.sizes[2,7:9])))
  i.tab[[9]] <- test.BLACK <- prop.test(tab.int.sizes[,10],c(sum(tab.int.sizes[1,10:11]),sum(tab.int.sizes[2,10:11])))
  i.tab[[10]] <- test.NotBLACK <- prop.test(tab.int.sizes[,11],c(sum(tab.int.sizes[1,10:11]),sum(tab.int.sizes[2,10:11])))
  i.tab[[11]] <- test.Chronic <- prop.test(tab.int.sizes[,12],c(sum(tab.int.sizes[1,12:13]),sum(tab.int.sizes[2,12:13]))) #sig
  i.tab[[12]] <- test.NotChronic <- prop.test(tab.int.sizes[,13],c(sum(tab.int.sizes[1,12:13]),sum(tab.int.sizes[2,12:13]))) #sig
  i.tab[[13]] <- test.FLUA <- prop.test(tab.int.sizes[,14],
                         c(sum(!is.na(hhls.in.prim.3$FluType) & hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & hhls.in.prim.3$Intervention..y.n. == "No",na.rm=T),
                           sum(!is.na(hhls.in.prim.3$FluType) & hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & hhls.in.prim.3$Intervention..y.n. == "Yes",na.rm=T))) #close
  age.cont <- hhls.in.prim.3[which(hhls.in.prim.3$Intervention..y.n. == "No" & hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & !is.na(hhls.in.prim.3$Relationship.to.Sick.Child) & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,"Age.of.Family.Member"]
  age.int <- hhls.in.prim.3[which(hhls.in.prim.3$Intervention..y.n. == "Yes" & hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & !is.na(hhls.in.prim.3$Relationship.to.Sick.Child) & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,"Age.of.Family.Member"]
  test.MainChildAge <- t.test(log(age.cont),log(age.int))
  i.tab[[14]] <- wilcox.test(age.cont,age.int)

  test.MeanHHLSize <- t.test(num.cont,num.int)
  i.tab[[15]] <- wilcox.test(num.cont,num.int)
  #which have a p-value less than 0.05
  p.less.than.05.int <- names(tab.int.sizes)[c(which(lapply(i.tab,function(x) x$p.value) < .05)+1)] #shift by one to account for first column
  print(paste("p < 0.05 (by intervention):",p.less.than.05.int))

  tab.type <- ddply(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids,],.(FluType),summarise,
                    "NoFluShot"=mean(Immunization.of.Family.Member == "No flu shot",na.rm=T),
                    "FluShotThisYear"=mean(Immunization.of.Family.Member == "This school year",na.rm=T),
                    "FluShotLastYear"=mean(Immunization.of.Family.Member == "Last school year",na.rm=T),
                    "PercentU15"=mean(Age.of.Family.Member <= 15,na.rm=T),
                    "PercentO15"=mean(Age.of.Family.Member > 15,na.rm=T),
                    "WorkFT"=mean(Activity.Out.of.Home.Fam.Member=="Work FT",na.rm=T),
                    "InSchoolorDayCare"=mean(Activity.Out.of.Home.Fam.Member=="School" | Activity.Out.of.Home.Fam.Member=="Daycare",na.rm=T),
                    "NotSchoolDayCareWorkFT"=mean(Activity.Out.of.Home.Fam.Member !="School" & Activity.Out.of.Home.Fam.Member !="Daycare",na.rm=T),
                    "Black"=mean(Race.of.Family.Member == "Black",na.rm=T),
                    "NotBlack"=mean(Race.of.Family.Member != "Black",na.rm=T),
                    "Chronic"=mean(Chronic.Illness..Fam.Member == "Yes",na.rm=T),
                    "NotChronic"=mean(Chronic.Illness..Fam.Member != "Yes",na.rm=T),
                    "MainChildAge"=mean(Age.of.Family.Member[Relationship.to.Sick.Child == "Sick Child (Self)"],na.rm=T))


  tab.type.sizes <- ddply(hhls.in.prim.3[hhls.in.prim.3$Sick.Childs.ID %in% ids.4,],.(FluType),summarise,
                          "NoFluShot"=sum(Immunization.of.Family.Member == "No flu shot",na.rm=T),
                          "FluShotThisYear"=sum(Immunization.of.Family.Member == "This school year",na.rm=T),
                          "FluShotLastYear"=sum(Immunization.of.Family.Member == "Last school year",na.rm=T),
                          "PercentU15"=sum(Age.of.Family.Member <= 15,na.rm=T),
                          "PercentO15"=sum(Age.of.Family.Member > 15,na.rm=T),
                          "WorkFT"=sum(Activity.Out.of.Home.Fam.Member=="Work FT",na.rm=T),
                          "InSchoolorDayCare"=sum(Activity.Out.of.Home.Fam.Member=="School" | Activity.Out.of.Home.Fam.Member=="Daycare",na.rm=T),
                          "NotSchoolDayCareWorkFT"=sum(Activity.Out.of.Home.Fam.Member !="School" & Activity.Out.of.Home.Fam.Member !="Daycare",na.rm=T),
                          "Black"=sum(Race.of.Family.Member == "Black",na.rm=T),
                          "NotBlack"=sum(Race.of.Family.Member != "Black",na.rm=T),
                          "Chronic"=sum(Chronic.Illness..Fam.Member == "Yes",na.rm=T),
                          "NotChronic"=sum(Chronic.Illness..Fam.Member != "Yes",na.rm=T),
                          "MainChildAge"=sum(Relationship.to.Sick.Child == "Sick Child (Self)" & !is.na(Age.of.Family.Member) ,na.rm=T))

  num.A <- ddply(hhls.in.prim.3[which(hhls.in.prim.3$FluType == "A" & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,],.(Sick.Childs.ID),function(x) x$Number.of.people.in.family[1])[,2]
  num.B <- ddply(hhls.in.prim.3[which(hhls.in.prim.3$FluType == "B" & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,],.(Sick.Childs.ID),function(x) x$Number.of.people.in.family[1])[,2]
  tab.type <- cbind(tab.type,"Mean Household Size" = c(mean(num.A),mean(num.B)))

                                        #prop tests
  t.tab <- vector("list",14)

  t.tab[[1]] <- test.NOFLUSHOT.type <- prop.test(tab.type.sizes[,2],c(sum(tab.type.sizes[1,2:4]),sum(tab.type.sizes[2,2:4])))
  t.tab[[2]] <- test.FluShotThisYear.type <-prop.test(tab.type.sizes[,3],c(sum(tab.type.sizes[1,2:4]),sum(tab.type.sizes[2,2:4])))
  t.tab[[3]] <- test.FluShotLastYear.type <-prop.test(tab.type.sizes[,4],c(sum(tab.type.sizes[1,2:4]),sum(tab.type.sizes[2,2:4])))
  t.tab[[4]] <- test.U15.type <- prop.test(tab.type.sizes[,5],c(sum(tab.type.sizes[1,5:6]),sum(tab.type.sizes[2,5:6])))
  t.tab[[5]] <- test.O15.type <- prop.test(tab.type.sizes[,6],c(sum(tab.type.sizes[1,5:6]),sum(tab.type.sizes[2,5:6])))
  t.tab[[6]] <- test.WorkFT.type <- prop.test(tab.type.sizes[,7],c(sum(tab.type.sizes[1,7:9]),sum(tab.type.sizes[2,7:9])))
  t.tab[[7]] <- test.INSCHOOLorDAYCARE.type <- prop.test(tab.type.sizes[,8],c(sum(tab.type.sizes[1,7:9]),sum(tab.type.sizes[2,7:9])))
  t.tab[[8]] <- test.NotSchoolDayCareWorkFT.type <- prop.test(tab.type.sizes[,9],c(sum(tab.type.sizes[1,7:9]),sum(tab.type.sizes[2,7:9])))
  t.tab[[9]] <- test.BLACK.type <- prop.test(tab.type.sizes[,10],c(sum(tab.type.sizes[1,10:11]),sum(tab.type.sizes[2,10:11])))
  t.tab[[10]] <- test.NotBLACK.type <- prop.test(tab.type.sizes[,11],c(sum(tab.type.sizes[1,10:11]),sum(tab.type.sizes[2,10:11])))
  t.tab[[11]] <- test.Chronic.type <- prop.test(tab.type.sizes[,12],c(sum(tab.type.sizes[1,12:13]),sum(tab.type.sizes[2,12:12])))
  t.tab[[12]] <- test.NotChronic.type <- prop.test(tab.type.sizes[,13],c(sum(tab.type.sizes[1,12:13]),sum(tab.type.sizes[2,12:13])))

  age.A <- hhls.in.prim.3[which(hhls.in.prim.3$FluType == "A" & hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & !is.na(hhls.in.prim.3$Relationship.to.Sick.Child) & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,"Age.of.Family.Member"]
  age.B <- hhls.in.prim.3[which(hhls.in.prim.3$FluType == "B" & hhls.in.prim.3$Relationship.to.Sick.Child == "Sick Child (Self)" & !is.na(hhls.in.prim.3$Relationship.to.Sick.Child) & hhls.in.prim.3$Sick.Childs.ID %in% ids) ,"Age.of.Family.Member"]
  test.MainChildAge <- t.test(log(age.A),log(age.B))
  t.tab[[13]] <- wilcox.test(age.A,age.B)
  test.MeanHHLSize <- t.test(num.A,num.B)
  t.tab[[14]] <-wilcox.test(num.A,num.B)
  p.less.than.05.type <- names(tab.type.sizes)[c(which(lapply(t.tab,function(x) x$p.value) < .05)+1)] #shift by one to account for first column
  print(paste("p < 0.05 (by type):",p.less.than.05.type))
  return(list("tab.int"=tab.int,"tab.int.sizes"=tab.int.sizes,"i.tab"=i.tab,"tab.type"=tab.type,"tab.type.sizes"=tab.type.sizes,"t.tab"=i.tab))
}

#make new dataframe for analysis
make.gee.df <- function(secs.gee){
gee.df <- matrix(nrow=nrow(secs.gee),ncol=16)
colnames(gee.df) <- c("ID","Age","HouseholdSize","Gender0f1m","Black0n1y","Chronic0n1y","SchoolOrDaycare0n1y","OnlyResid0n1y","ImmunThisYear0n1y","Cont0Int1","Grade","Flu0A1B","SICK.NoTemp","SICK.Temp","Grade.v2","Age.cat")
gee.df <- as.data.frame(gee.df)
gee.df$ID <- secs.gee$Sick.Childs.ID
gee.df$Age <- secs.gee$Age.of.Family.Member
gee.df$HouseholdSize <- secs.gee$Number.of.people.in.family
gee.df$Gender0f1m <- ifelse(secs.gee$Gender.of.Family.Member == "Female",0,1)
gee.df$Black0n1y <- ifelse(secs.gee$Race.of.Family.Member == "Black",1,0)
gee.df$Chronic0n1y <- ifelse(secs.gee$Chronic.Illness..Fam.Member == "No",0,1)
gee.df$SchoolOrDaycare0n1y <- ifelse(secs.gee$Activity.Out.of.Home.Fam.Member == "School" | secs.gee$Activity.Out.of.Home.Fam.Member == "Daycare",1,0)
gee.df$OnlyResid0n1y <- ifelse(secs.gee$Only.Residence.of.Fam.Member == "No",0,1)
gee.df$ImmunThisYear0n1y <- ifelse(secs.gee$Immunization.of.Family.Member == "This school year",1,0)
gee.df$Cont0Int1 <- ifelse(secs.gee$Intervention..y.n. == "No",0,1)
gee.df$Grade <- secs.gee$Grade

gee.df$Flu0A1B <- ifelse(secs.gee$FluType == "A",0,1)
gee.df$SICK.NoTemp <- secs.gee$SICK.NoTemp
gee.df$SICK.Temp <- secs.gee$SICK.Temp
for (i in 1:nrow(gee.df)){
  if (gee.df$Grade[i] == 0 | gee.df$Grade[i] == 1){
    gee.df$Grade.v2[i] <- 0
  } else if (gee.df$Grade[i] == 2 | gee.df$Grade[i] == 3){
    gee.df$Grade.v2[i] <- 1
  } else if (gee.df$Grade[i] == 4 | gee.df$Grade[i] == 5){
    gee.df$Grade.v2[i] <- 2
  }
}
for (i in 1:nrow(gee.df)){
  if(is.na(gee.df$Age[i])){
    gee.df$Age.cat[i] <- NA
  } else if(gee.df$Age[i] < 5){
    gee.df$Age.cat[i] <- 1
  } else if(gee.df$Age[i] >=5 & gee.df$Age[i] < 19) {
    gee.df$Age.cat[i] <- 2
  } else if(gee.df$Age[i] >=19 & gee.df$Age[i] <=50){
    gee.df$Age.cat[i] <- 3
  } else if(gee.df$Age[i] >= 51){
    gee.df$Age.cat[i] <- 4
  }
}
  gee.df$Age.cat <- factor(gee.df$Age.cat,levels=c(3,1,2,4)) #make the middle cat the reference group
  return(gee.df)
}
