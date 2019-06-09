geteval = function(Delta_real, Delta_est, p){
    diag(Delta_real) = 0; diag(Delta_est) = 0
    Deltabias = norm(Delta_real - Delta_est, type = "F")
    
    ind = upper.tri(Delta_real,  diag = FALSE)
    TPR = sum((Delta_real != 0 & Delta_est != 0)[ind]) / sum(Delta_real[ind] != 0)
    TNR = sum((Delta_real == 0 & Delta_est == 0)[ind]) / sum(Delta_real[ind] == 0)
    
    return(data.frame(Deltabias, TPR, TNR))
}

difpathadd = function(x){
  len = length(x)
  p = nrow(x[[1]])
  difidx = lapply(x, function(mat){
      ret = (mat!=0)*1
      diag(ret) = 0
      return(ret)
  })
  difidx[[len+1]] = matrix(0,p,p)
  ret = matrix(1,p,p); diag(ret) = 0
  difidx[[len+2]] = ret
  return(difidx)
}


getdifpath = function(path, path_new){
    path = lapply(path, function(x)return((x+t(x))/2))
    path_new = lapply(path_new, function(x)return((x+t(x))/2))
    len = length(path)
    difidx = lapply(1:len, function(k){
        ret = (as.matrix(path[[k]]) - as.matrix(path_new[[k]])!=0)*1
        diag(ret) = 0
        return(ret)
    })
    
    p = nrow(path[[1]])
    difidx[[len+1]] = matrix(0,p,p)
    ret = matrix(1,p,p); diag(ret) = 0
    difidx[[len+2]] = ret
    return(difidx)
}

onesim = function(Sigma, Theta, mu, p, n, seed, model){
  set.seed(seed)
  data = exp(mvrnorm(n = n, mu = mu, Sigma = Sigma))
  data = data / rowSums(data)
  clr.data = log(data) %*% (diag(p) - 1/p * matrix(1, p, p))
  
  data_new = exp(mvrnorm(n = n, mu = mu_new, Sigma = Sigma_new))
  data_new = data_new / rowSums(data_new)
  clr.data_new = log(data_new) %*% (diag(p) - 1/p * matrix(1, p, p))
  
  Delta = Theta- Theta_new
  Delta_offdiag = 1*(Delta != 0)
  diag(Delta_offdiag) = 0
  
  ##se(mb)
  mb.stime = Sys.time()
  mb.res = sparseiCov(clr.data, method = 'mb', lambda.min.ratio = 0.0001, nlambda = 50, sym = "or")
  mb.res_new = sparseiCov(clr.data_new, method = 'mb', lambda.min.ratio = 0.0001, nlambda = 50, sym = "or")
  mb.roc = huge.roc(getdifpath(mb.res$beta, mb.res_new$beta), Delta_offdiag, verbose = FALSE)
  mb.roc

  mb.res = sparseiCov(clr.data, method = 'mb', lambda.min.ratio = 0.1, nlambda = 20, sym = "or")
  mb.res_new = sparseiCov(clr.data_new, method = 'mb', lambda.min.ratio = 0.1, nlambda = 20, sym = "or")
  mb.res = icov.select(mb.res, criterion = 'stars')
  mb.res_new = icov.select(mb.res_new, criterion = 'stars')
  mb.icov = symBeta(getOptBeta(mb.res), mode='ave')
  mb.icov_new = symBeta(getOptBeta(mb.res_new), mode='ave')
  mbval = geteval(Delta_real = Delta, Delta_est = as.matrix(mb.icov) - as.matrix(mb.icov_new), p = p)
  mb.etime = Sys.time()

  ##se(gl)
  gl.stime = Sys.time()
  gl.res = sparseiCov(clr.data, method = 'glasso', lambda.min.ratio = 0.0001, nlambda = 50)
  gl.res_new = sparseiCov(clr.data_new, method = 'glasso', lambda.min.ratio = 0.0001, nlambda = 50)
  gl.roc = huge.roc(getdifpath(gl.res$icov, gl.res_new$icov), Delta_offdiag, verbose = FALSE)
  gl.roc

  gl.res = sparseiCov(clr.data, method = 'glasso', lambda.min.ratio = 0.001, nlambda = 15)
  gl.res_new = sparseiCov(clr.data_new, method = 'glasso', lambda.min.ratio = 0.001, nlambda = 15)
  gl.res = icov.select(gl.res, criterion = 'stars')
  gl.res_new = icov.select(gl.res_new, criterion = 'stars')
  gl.icov = gl.res$opt.icov; gl.icov_new = gl.res_new$opt.icov
  glval = geteval(Delta_real = Delta, Delta_est = as.matrix(gl.icov) - as.matrix(gl.icov_new), p = p)
  gl.etime = Sys.time()


  ###result of fang
  fang.stime = Sys.time()
  fang = gcoda(x = data, counts = F, lambda.min.ratio = 1e-3,
               nlambda = 50, ebic.gamma = 0.5)
  fang_new = gcoda(x = data_new, counts = F, lambda.min.ratio = 1e-3,
               nlambda = 50, ebic.gamma = 0.5)
  fang.roc = huge.roc(getdifpath(fang$icov, fang_new$icov), Delta_offdiag, verbose = FALSE)
  fang.roc
  fangval = geteval(Delta_real = Delta, Delta_est = as.matrix(fang$opt.icov) - as.matrix(fang_new$opt.icov), p = p)
  fang.etime = Sys.time()


  ##CDtrace approximate
  cdtrace.stime_app = Sys.time()
  cdtrace.res_app = cdtrace_path(data, exact = FALSE, rho = 1000, eps = 1e-8, lambda.min.ratio = 1e-2, nlambda = 50,
                             tol = 1e-3, itermax = 1000, fastitermax = 1000, fasttol = 1e-3)
  cdtrace.res_app_new = cdtrace_path(data_new, exact = FALSE, rho = 1000, eps = 1e-8, lambda.min.ratio = 1e-2, nlambda = 50,
                                 tol = 1e-3, itermax = 1000, fastitermax = 1000, fasttol = 1e-3)
  cdtrace.roc_app = huge.roc(getdifpath(cdtrace.res_app$icovpath, cdtrace.res_app_new$icovpath), Delta_offdiag, verbose = FALSE)
  cdtrace.roc_app
  cdtraceval_app = geteval(Delta_real = Delta, Delta_est = cdtrace.res_app$icov - cdtrace.res_app_new$icov, p = p)
  cdtrace.etime_app = Sys.time()

  ##CDtrace
  cdtrace.stime = Sys.time()
  cdtrace.res = cdtrace_path(data, exact = TRUE, rho = 1000, eps = 1e-8, lambda.min.ratio = 1e-2, nlambda = 50,
                          tol = 1e-3, itermax = 1000, fastitermax = 1000, fasttol = 1e-3)
  cdtrace.res_new = cdtrace_path(data_new, exact = TRUE, rho = 1000, eps = 1e-8, lambda.min.ratio = 1e-2, nlambda = 50,
                             tol = 1e-3, itermax = 1000, fastitermax = 1000, fasttol = 1e-3)
  cdtrace.roc = huge.roc(getdifpath(cdtrace.res$icovpath, cdtrace.res_new$icovpath), Delta_offdiag, verbose = FALSE)
  cdtrace.roc
  cdtraceval = geteval(Delta_real = Delta, Delta_est = cdtrace.res$icov - cdtrace.res_new$icov, p = p)
  cdtrace.etime = Sys.time()

  
  ##DPM
  dpm.stime = Sys.time()
  dpm.res = dpm(clr.data, clr.data_new, lambda=NULL, nlambda=50,lambda.min.ratio=0.01,
                 rho=NULL,shrink=NULL,prec=0.0001,max.ite=5,
                 correlation=FALSE,perturb=FALSE,tuning=c("bic"))
  dpm.roc = huge.roc(difpathadd(dpm.res$dpm), Delta_offdiag, verbose = FALSE)
  dpm.roc
  dpmval = geteval(Delta_real = Delta, Delta_est = dpm.res$dpm[[dpm.res$opt[1]]], p = p)
  dpmval
  dpm.etime = Sys.time()
  
  ##yuan
  yuan.stime = Sys.time()
  yuan = Dpmdtl(clr.data, clr.data_new, lambda=NULL, nlambda=50, lambda.min.ratio=0.001,
                    rho=NULL,shrink=NULL,prec=0.001, correlation=FALSE,tuning=c("bic"))
  yuan.roc = huge.roc(difpathadd(yuan$Dpmdtl), Delta_offdiag, verbose = FALSE)
  yuan.roc
  yuanval = geteval(Delta_real = Delta, Delta_est = yuan$Dpmdtl[[yuan$opt[1]]], p = p)
  yuanval
  yuan.etime = Sys.time()
  
  
  ##FGL
  fgl.stime = Sys.time()
  fgl.res = lapply(yuan$lambda, function(ll){
      ret = JGL(Y = list(clr.data, clr.data_new),penalty="fused",return.whole.theta=TRUE,
                lambda1 = ll ,lambda2 = 0.01)
      return(ret$theta[[1]] - ret$theta[[2]])
  })
  fgl.roc = huge.roc(difpathadd(fgl.res), Delta_offdiag, verbose = FALSE)
  fgl.roc
  fglval = geteval(Delta_real = Delta, Delta_est = fgl.res[[yuan$opt[1]]], p = p)
  fglval
  fgl.etime = Sys.time()
  
  ##GGL
  ggl.stime = Sys.time()
  ggl.res = lapply(yuan$lambda, function(ll){
      ret = JGL(Y = list(clr.data, clr.data_new),penalty="group",return.whole.theta=TRUE,
                lambda1 = ll ,lambda2 = 0.01)
      return(ret$theta[[1]] - ret$theta[[2]])
  })
  ggl.roc = huge.roc(difpathadd(ggl.res), Delta_offdiag, verbose = FALSE)
  ggl.roc
  gglval = geteval(Delta_real = Delta, Delta_est = ggl.res[[yuan$opt[1]]], p = p)
  gglval
  ggl.etime = Sys.time()
  
  
  maxF1_score = data.frame(mb = max(mb.roc$F1), gl = max(gl.roc$F1), 
                           fang = max(fang.roc$F1), cdtrace_app = max(cdtrace.roc_app$F1),
                           cdtrace = max(cdtrace.roc$F1), dpm = max(dpm.roc$F1),
                           yuan = max(yuan.roc$F1),fgl = max(fgl.roc$F1),ggl = max(ggl.roc$F1))
  AUCtab = data.frame(mb = mb.roc$AUC, gl = gl.roc$AUC, 
                      fang = fang.roc$AUC, cdtrace_app = cdtrace.roc_app$AUC, 
                      cdtrace = cdtrace.roc$AUC, dpm = dpm.roc$AUC,
                      yuan = yuan.roc$AUC, fgl = fgl.roc$AUC, ggl = ggl.roc$AUC)
  valtab = rbind(mbval, glval, fangval, cdtraceval_app, cdtraceval,dpmval,yuanval,fglval,gglval)
  rownames(valtab) = c("mb", "gl", "fang", "cdtrace_app", "cdtrace","dpm","yuan","fgl","ggl")
  timetab = data.frame(mb = mb.etime - mb.stime,
                       gl = gl.etime - gl.stime,
                       fang = fang.etime - fang.stime,
                       cdtrace_app = cdtrace.etime_app - cdtrace.stime_app,
                       cdtrace = cdtrace.etime - cdtrace.stime,
                       dpm = dpm.etime - dpm.stime,
                       yuan = yuan.etime - yuan.stime,
                       fgl = fgl.etime - fgl.stime,
                       ggl = ggl.etime - ggl.stime)
  
  save(seed, p, n, valtab, maxF1_score, AUCtab, timetab,
       fang.roc, cdtrace.roc, mb.roc, gl.roc, cdtrace.roc_app,dpm.roc,yuan.roc,fgl.roc,ggl.roc,
       file = paste0(model, "-p", p, "-n", n, "-seed", seed, ".rdata"))
  return(NULL)
}
