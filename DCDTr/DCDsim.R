######assumption not right####
starttime = Sys.time()

path = "/public/home/heshun/hs/DCDtrace"
setwd(path)
source("gCoda-master//R//gcoda.R")
require(huge)
source("simfunc.R")

require(doSNOW)
require(parallel)
cl = makeSOCKcluster(rep("localhost", times = 10))
registerDoSNOW(cl)

p = 50
n = 400
#model = "random"
#model = "neighbor"
#model = "block"
#model = "hub" 
#model = "scale-free"
#model = "band"
model = "cluster"
for(model in c("random", "hub", "neighbor", "block", "band", "scale-free")){#
    if(model == "random"){
        g = NULL; prob = 0.5; upper = 0.4; lower = 0.2
    }else if(model == "hub"){
        g = 3; prob = NULL; upper = 0.4; lower = 0.2
    }else if(model == "cluster"){
        g = 3; prob = 0.6; upper = 0.4; lower = 0.2
    }else if(model == "band"){
        g = 4; prob = NULL; upper = NULL; lower = NULL
    }else if(model == "scale-free"){
        g = NULL; prob = NULL; upper = 0.4; lower = 0.2
    }
    if(model == "neighbor"){
        set.seed(1)
        Theta = matrix(0,p,p)
        locmat = matrix(0, p, p)
        d  = runif(n = p, 0, 1)
        for(i in 1:p){
            x = abs(d[i] - d)
            loc = order(x)[2:11]
            locmat[i, loc] = 1; locmat[loc, i] = 1
        }
        for(i in 1:(p-1)){
            for(j in (i+1):p){
                if(locmat[i,j] == 1){
                    val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
                        runif(1, 0.2, 0.4)
                    Theta[i, j] = val; Theta[j, i] = val
                }
            }
        }
    }else if(model == "block"){
        set.seed(1)
        Theta = matrix(0,p,p)
        block = sample(p, size = p, replace = FALSE)
        blocksize = p / 5
        block = lapply(1:5, function(i){
            sort(block[(blocksize*(i-1)+1):(blocksize*i)])
        })
        
        for(i in 1:5){
            for(j in block[[i]]){
                for(k in block[[i]][block[[i]] > j]){
                    if(rbinom(n = 1, size = 1, prob = 0.5)){
                        val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
                            runif(n = 1, 0.2, 0.4)
                        Theta[k,j] = val
                        Theta[j,k] = val
                    }}}}
        
        for(i in 1:4){
            for(j in block[[i]]){
                for(k in unlist(lapply(i:5, function(l)block[[l]]))){
                    if(rbinom(n = 1, size = 1, prob = 0.3)){
                        val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
                            runif(n = 1, 0.2, 0.4)
                        Theta[k,j] = val
                        Theta[j,k] = val
                    }}}}
    }else{
        ##generate network
        set.seed(1)
        net_generator = huge.generator(n = 2*p, d = p, graph = model, v = NULL, u = NULL, 
                                       g = g, prob = prob, vis = FALSE, verbose = FALSE)
        #plot(net_generator, align = F)
        net = as.matrix(net_generator$theta)
        
        if(model == "band"){
            link_strength = matrix(0,p,p)
            for(i in 1:(p-1))link_strength[i,i+1] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.2, 0.25)
            for(i in 1:(p-2))link_strength[i,i+2] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.15, 0.2)
            for(i in 1:(p-3))link_strength[i,i+3] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.1, 0.15)
            for(i in 1:(p-4))link_strength[i,i+4] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.05, 0.1)
        }else{
            link_strength = matrix((2 * rbinom(n = p^2, size = 1, prob = 0.5) - 1) * 
                                       runif(n = p^2, lower, upper), p, p)
        }
        link_strength[lower.tri(link_strength, diag = TRUE)] =  0
        link_strength = link_strength + t(link_strength)
        
        Theta = net * link_strength
    }
    
    nzeroloc = NULL
    zeroloc = NULL
    for(i in 1:(p-1)){
        for(j in (i+1):p){
            if(Theta[i,j]!=0){
                nzeroloc = rbind(nzeroloc, c(i,j))
            }else{
                zeroloc = rbind(zeroloc, c(i,j))
            }
        }  
    }
    
    if(model %in% c("random", "neighbor", "block","band")){
        change_ratio = 0.1
    }else if(model %in% c("hub" ,"scale-free")){
        change_ratio = 0.4
    }else if(model %in% c("cluster")){
        change_ratio = 0.2
    }
    changenum = round(nrow(nzeroloc)*change_ratio,0)
    changeidx1 = sample(1:nrow(nzeroloc), changenum, replace=FALSE)
    changeidx2 = sample(1:nrow(zeroloc), changenum, replace=FALSE)
    Theta_new = Theta
    for(i in 1:changenum){
        ii = nzeroloc[changeidx1[i],1]
        jj = nzeroloc[changeidx1[i],2]
        Theta_new[ii, jj] = Theta_new[jj,ii] = 0
        ii = zeroloc[changeidx2[i],1]
        jj = zeroloc[changeidx2[i],2]
        Theta_new[ii, jj] = Theta_new[jj,ii] = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.2, 0.4)
    }
    
    toadd = max(abs(min(eigen(Theta)$values)), abs(min(eigen(Theta_new)$values)))
    diag(Theta) = toadd + 0.3
    diag(Theta_new) = toadd + 0.3
    
    Sigma = solve(Theta); Sigma_new = solve(Theta_new)
    mu = runif(n = p, min  = -0.5, max = 0.5)
    mu_new = runif(n = p, min  = -0.5, max = 0.5)
    
    GG = diag(p) - 1/p*matrix(1,p,p)
    asbias = c(norm(GG %*% Theta - Theta %*% GG, 'F'), norm(GG %*% Sigma - Sigma %*% GG, 'F'))
    asbias_new = c(norm(GG %*% Theta_new - Theta_new %*% GG, 'F'), norm(GG %*% Sigma_new - Sigma_new %*% GG, 'F'))
    asbias; asbias_new
    
    setwd(path)
    modelpath = paste0(path, "//", model)
    dir.create(modelpath, showWarnings = FALSE)
    setwd(modelpath)
    
    save(Theta, Sigma, mu, asbias,Theta_new, Sigma_new, mu_new, asbias_new,
         file = paste0(model, "-p", p, "-net.rdata"))
    
    
    for(n in c(100, 200, 300, 400)){# 
        cat("graph:", model, "  n:", n, "\n")
        void = foreach(seed = 1:120)%dopar%{
            path1 = paste0(path, "/SpiecEasi-master/R")
            setwd(path1)
            for(f in list.files(path1)){
                source(f)
            }
            path2 = paste0(path, "/CDtrace")
            setwd(path2)
            for(f in list.files(path2)){
                source(f)
            }
            setwd(path)
            source("gCoda-master//R//gcoda.R")
            require(huge)
            source("simfunc.R")
            require(Difdtl)
            #require(codaJGL)
            require(JGL)
            path3 = paste0(path, "/dpm-master")
            setwd(path3)
            source("dpm.R")
            dyn.load("dpm.so")
            
            setwd(modelpath)
            errind = try(onesim(Sigma = Sigma, Theta = Theta, mu = mu, 
                                p = p, n = n, seed = seed, model = model),
                         silent = TRUE)
            
            if(inherits(errind, "try-error")){
                return(seed)
            }else{
                return(0)
            }
        }
        
        setwd(modelpath)
        save(void, file = paste0(model, "-p", p, "-n", n, "-errseed.rdata"))
    }
}

stopCluster(cl)
endtime = Sys.time()

