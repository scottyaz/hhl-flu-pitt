## Effect Size Simulations for Household Transmission of Influenza A and B in a School-based Study of Non-Pharmaceutical Interventions published in the Journal Epidemics
## A. Azman et al 2013, Epidemics

## Function to geneate final attack sizes of a number of households and then estimate the MLEs
run.effect.size.sims <- function(q.cont.vec=rep(.75,20),
                                 q.int.vec=seq(.85,.97,length=20),
                                 alpha,n.hhls=89,
                                 num.conts=44,
                                 n.sims=1000){

    sim.res <- array(dim=c(n.sims,4,length(q.cont.vec))) #for final results
    ##get household sizes and emperical dsitribution
    household.sizes <- table(ddply(hhls.in.prim.3,.(Sick.Childs.ID),function(x) mean(x$Number.of.people.in.family))[,2])
    emp.hhl.dist <- household.sizes/sum(household.sizes)

    ## this final outbreak size probabilities and household sizes and gives a final outbreak size for the sample
    get.final.tab <- function(sim.row,ps){
        sim.tab <- table(factor(sim.row,levels=2:6))  ##table of household sizes
        hhls.tab.temp <- matrix(nrow=6,ncol=5)
        for (k in 2:6) {
            if(k == 6){
                hhls.tab.temp[,k-1] <- table(factor(sample(0:(k-1),size=sim.tab[paste(k)],prob=ps[1:k,k-1],replace=T),levels=0:(k-1)))
            } else {
                add.zeros <- rep(0,6-k)
                tab.temp <- table(factor(sample(0:(k-1),sim.tab[paste(k)],prob=ps[1:k,k-1],replace=T),0:(k-1)))
                hhls.tab.temp[,k-1] <- c(tab.temp,add.zeros)
            }
        }

        colnames(hhls.tab.temp) <- 2:6
        rownames(hhls.tab.temp) <- 0:5
        return(hhls.tab.temp)
    }

    ## get probabilities for each household size based on truncated distribution
    make.prob.matrix <- function(q,alpha){
        ps <- NULL
        for (i in 2:6){
            p.temp <- fs_probs_mod3(n=i,
                                    q=q,
                                    alpha=alpha,
                                    total.vector.length = 7,
                                    log=FALSE,
                                    truncated=FALSE)
            ps <- cbind(ps,p.temp)
        }
        return(ps)
    }

    for (j in 1:length(q.cont.vec)){
        q.cont <- q.cont.vec[j]
        q.int <- q.int.vec[j]
        ## matrix of hosuehold sizes for each simulation
        sim.hhl.sizes <- matrix(nrow=n.sims,ncol=n.hhls)
        int.assignments.sims <- matrix(nrow=n.sims,ncol=n.hhls)
        for (i in 1:n.sims) {
            sim.hhl.sizes[i,] <- sample(2:6,n.hhls,emp.hhl.dist,replace=T)
            int.assignments.sims[i,] <- sample(c(0,1),n.hhls,prob=c(.5,.5),replace=T)
        }

        ps.cont <- make.prob.matrix(q.cont,alpha)
        ps.int <- make.prob.matrix(q.int,alpha)
                                        #go through each household and decide who gets sick based on final size probs
                                        #sim.row - vector of household sizes
                                        #split into control and intervention households
        conts <- sim.hhl.sizes[,1:num.conts]
        ints <- sim.hhl.sizes[,(num.conts+1):ncol(sim.hhl.sizes)]

        ##generate final outbreak size for control and intervention groups
        cat(sprintf("Generating final outbreak sizes for %f and %f \n",q.cont,q.int))
        tables.cont = array(dim=c(6,5,n.sims))
        tables.int = array(dim=c(6,5,n.sims))
        for (i in 1:n.sims){
            tables.cont[,,i] <- get.final.tab(conts[i,],ps.cont)
            tables.int[,,i] <- get.final.tab(conts[i,],ps.int)
        }
                                        #get mles and cis for each simulation
                                        #sim.results <- matrix(nrow=n.sims,ncol=6)
        sim.results <- matrix(nrow=n.sims,ncol=4)
        colnames(sim.results) <- c("Q.cont","Q.int","LRTs","Sig")
        print("Estimating MLEs and Profile Likelihoods \n")
        for (i in 1:n.sims){
            if (i %% 10 == 0) cat(".")
            hhltab1 <- tables.cont[,,i]
            hhltab2 <- tables.int[,,i]
            flu.tables <- list(hhltab1,hhltab2)
            tmp <- likelihood.ratio.tests2(flu.tables)
            sim.results[i,] <- tmp
        }
        sim.res[,,j] <- sim.results
    }
    sim.overview <- apply(sim.res,3,function(x) mean(x[,4]))
                                        # return(sim.res)
    return(cbind(q.cont.vec,q.int.vec,sim.overview))
}

##probabiliies of escape for generating final outbreaks
q.cont.vec.75 <- rep(4^(-.71)*qlogis(.75),15)    #for control group
q.int.vec.75 <- seq(4^(-.71)*qlogis(.80),4^(-.71)*qlogis(.90),length=15)    #for intervention group
q.cont.vec.80 <- rep(4^(-.71)*qlogis(.80),15)
q.int.vec.80 <- seq(4^(-.71)*qlogis(.88),4^(-.71)*qlogis(.95),length=15)

effect.size.sims.75 <- run.effect.size.sims(q.cont.vec.75,
                                            q.int.vec.75,
                                            n.hhls=89,
                                            num.conts=44,
                                            n.sims=1000,
                                            alpha=-0.71)

effect.size.sims.80 <- run.effect.size.sims(q.cont.vec.80,
                                            q.int.vec.80,
                                            n.hhls=89,
                                            num.conts=44,
                                            n.sims=1000,
                                            alpha=-0.71)
# save data
dput(list(effect.size.sims.75,effect.size.sims.80),"Data/effect_size_sims.rda")

## how much would we have to reduce our SITP by the intervention to detect a change by a LRT?
sitp.target.80 <- approx(x=effect.size.sims.80[,3],
                         y=effect.size.sims.80[,2],
                         xout=0.80)$y

(.20 - (1 - plogis(sitp.target.80/4^(-0.71))))/.20
                                        #[1] 0.4835

sitp.target.75 <- approx(x=effect.size.sims.75[,3],
                         y=effect.size.sims.75[,2],
                         xout=0.75)$y

(.25 - (1 - plogis(sitp.target.75/4^(-0.71))))/.25
                                        # [1] 0.4166
