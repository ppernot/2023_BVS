
figDir = '../Figs'
tabDir = '../Tabs'
library(ErrViewLib)
library(CHNOSZ)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

fractHetero = function(S) {
  fHetero = c()
  for(i in seq_along(S)) {
    compo = CHNOSZ::count.elements(S[i])
    tot   = sum(compo, na.rm = TRUE)
    CandH = sum(compo['C'], compo['H'], na.rm = TRUE)
    fHetero[i] = 1 - CandH / tot
  }
  fHetero
}
scaleByLimits <- function(X, U, s, lwlims) {
  Us = U
  for(i in 1:length(s)) {
    if(i == 1)
      xmin = -Inf
    else
      xmin = lwlims[i-1]

    if(i == length(s))
      xmax  = Inf
    else
      xmax = lwlims[i]

    sel = X >= xmin & X < xmax
    if(sum(sel) > 0)
      Us[sel] = U[sel] * s[i]
  }
  return(Us)
}
fScale <- function(pars, uE, scaleMod, nBin, intrv) {
  # Assuming uE is already sorted
  if(scaleMod == 'bin-wise') {
    # Bin-wise scaling of uE
    for(i in 1:nBin) {
      sel  = intrv$lwindx[i]:intrv$upindx[i]
      uE[sel] = uE[sel] * pars[i]
    }

  } else if(scaleMod == 'element-wise') {
    # Element-wise scaling of uE
    for(i in seq_along(uE))
      uE[i] = uE[i] * pars[i]

  } else {
    # Linear scaling
    for(i in seq_along(uE))
      uE[i] = pars[1] * uE[i] + pars[2]

  }

  return(uE)
}
scoreFun0 <- function(Z, nBin, intrv) {
  AE = c()
  for(i in 1:nBin) {
    sel   = intrv$lwindx[i]:intrv$upindx[i]
    AE[i] = abs(log(mean(Z[sel]^2)))
  }
  return(mean(AE))
}
scoreFun1 <- function(E, uE, nBin, intrv) {
  MV = MSE = c()
  for(i in 1:nBin) {
    sel    = intrv$lwindx[i]:intrv$upindx[i]
    MV[i]  = mean(uE[sel]^2)
    MSE[i] = mean(E[sel]^2)
  }
  return(
    list(
      ence = mean(abs(sqrt(MV) - sqrt(MSE)) / sqrt(MV)),
      uce  = mean(abs(MV - MSE))
    )
  )
}
nllFun = function(E, uE) {
  calib = 0.5*mean( E^2/uE^2  )
  sharp = 0.5*mean( log(uE^2) )
  nll   = calib + sharp + 0.5 * 1.837877
  return(
    list(
      calib = calib,
      sharp = sharp,
      nll   = nll
    )
  )
}
scoreFun = function(
  pars,
  E, uE, X1, X2,
  nBinOpt, intrvOpt, nBinStat, intrvStat,
  scoreComp ='1111',
  scaleMod = c('bin-wise','element-wise','linear')) {

  scaleMod = match.arg(scaleMod)

  # Apply scaling model to uE
  uEs = fScale(pars, uE, scaleMod, nBinOpt, intrvOpt)
  Z   = E /uEs
  
  score = 0

  # Mean calibration error or NLL
  if (substr(scoreComp, 1, 1) == "1")
    score = score + abs(log(mean(Z ^ 2)))
  else if (substr(scoreComp, 1, 1) == "2")
    score = score + nllFun(E, uEs)$nll

  # Mean consistency error
  if (substr(scoreComp, 2, 2) == "1") {
    io = order(uEs)
    score = score + scoreFun0(Z[io], nBinStat, intrvStat)
  }

  # Mean adaptivity error
  ## X1
  if (substr(scoreComp, 3, 3) == "1") {
    io = order(X1, uEs)
    score = score + scoreFun0(Z[io], nBinStat, intrvStat)
  }
  ## X2
  if (substr(scoreComp, 4, 4) == "1") {
    io = order(X2, uEs)
    score = score + scoreFun0(Z[io], nBinStat, intrvStat)
  }

  return(score)
}
globScore <- function(E, uE, X1, X2, nBin, intrv) {
  
  M = length(uE)
  io0 = order(uE)
  io1 = order(X1,1:M)
  io2 = order(X2,1:M)
  
  Z = E / uE
  Scal  = abs(log(mean(Z^2)))
  Scon  = scoreFun0(Z[io0], nBin, intrv)
  SX1   = scoreFun0(Z[io1], nBin, intrv)
  SX2   = scoreFun0(Z[io2], nBin, intrv)
  Stot  = Scal + Scon + SX1 + SX2
  li    = scoreFun1(E[io0],uE[io0], nBin, intrv)
  ence  = li$ence
  uce   = li$uce
  nl    = nllFun(E,uE)
  nll   = nl$nll
  calib = nl$calib
  sharp = nl$sharp
  
  return(
    list(
      Scal  = Scal, 
      Scon  = Scon, 
      SX1   = SX1, 
      SX2   = SX2,
      Stot  = Stot,
      ence  = ence, 
      uce   = uce, 
      nll   = nll, 
      calib = calib, 
      sharp = sharp
    )
  )
}
optScaleFun <- function(
  Etrain, uEtrain, X1train, X2train,
  scaleMod, scoreComp, scoreFun,
  nBinOpt, intrvOpt, nBinStat, intrvStat
) {

  if(scaleMod == 'bin-wise') {
    s0 = c()
    for(i in 1:nBinOpt) {
      sel  = intrvOpt$lwindx[i]:intrvOpt$upindx[i]
      s0[i] = sqrt(mean((Etrain[sel]/uEtrain[sel])^2))
    }
    lb =  0.1 * s0
    ub = 10.0 * s0
    algorithm = "NLOPT_GN_DIRECT" #"NLOPT_LN_BOBYQA"
    maxeval = 1e5
    start = scoreFun(s0, Etrain, uEtrain, X1train, X2train,
                     nBinOpt, intrvOpt, nBinStat, intrvStat,
                     scoreComp, scaleMod)

  } else if(scaleMod == 'element-wise') {
    # c = sqrt(mean((Etrain/uEtrain)^2))
    s0 = abs(Etrain/uEtrain)
    lb =  0.01 * s0
    ub = 10.0  * s0
    algorithm =  "NLOPT_LN_BOBYQA" #"NLOPT_GN_DIRECT"
    maxeval = 1e5

  } else {
    # linear
    c = sqrt(mean((Etrain/uEtrain)^2))
    s0 = c()
    s0[1] = 1
    s0[2] = 0
    lb = c(0,-1)
    ub = c(2,1)
    algorithm = "NLOPT_GN_AGS"
    maxeval = 1e5
  }

  # Optimize
  cat('\n\n Global optim for ',scoreComp,'\n')
  resOpt = nloptr::nloptr(
    x0 = s0,
    eval_f = scoreFun,
    opts =  list(
      algorithm   = algorithm,
      xtol_rel    = 1.0e-6,
      print_level = 0,
      maxeval     = maxeval),
    lb = lb,
    ub = ub,
    E = Etrain, uE = uEtrain, X1 = X2train, X2 = X2train,
    nBinOpt = nBinOpt, intrvOpt = intrvOpt,
    nBinStat = nBinStat, intrvStat = intrvStat,
    scoreComp = scoreComp, scaleMod = scaleMod)
  cat('Status : ', resOpt$status, ', Objective =',resOpt$objective,'\n\n')
  s1 = resOpt$solution
  if( resOpt$objective > start) # Global search did not improve s0
    s1 = s0

  # Refine solution
  cat('Local optim for ',scoreComp,'\n')
  status = 5
  nbEval = 0
  while(status == 5 & nbEval < maxeval) {
    resOpt = nloptr::nloptr(
      x0 = s1,
      eval_f = scoreFun,
      opts =  list(
        algorithm   = "NLOPT_LN_BOBYQA",
        xtol_rel    = 1.0e-6,
        print_level = 0,
        maxeval     = maxeval/5),
      lb = lb,
      ub = ub,
      E = Etrain, uE = uEtrain, X1 = X2train, X2 = X2train,
      nBinOpt = nBinOpt, intrvOpt = intrvOpt,
      nBinStat = nBinStat, intrvStat = intrvStat,
      scoreComp = scoreComp, scaleMod = scaleMod)
    status = resOpt$status
    nbEval = nbEval + maxeval / 5
    s1 = resOpt$solution
    cat('Status : ', status, ', Objective =',resOpt$objective,'\n')
  }

  s   = resOpt$solution
  end = resOpt$objective

  return(
    list(
      s0  = s0,
      s   = s,
      end = end
    )
  )
}
plotRes <- function(Etrain, uEtrain, uEtrains, X1train, X2train, nBin,
                    method = 'stud', plot = FALSE) {

  # Control plots, but mostly to generate ZMS validation stats
  if(plot) {
    par(mfrow = c(1,2))
    plot(uEtrain, Etrain); grid()
    abline(a=0, b=-2, lty=2, col=2)
    abline(a=0, b=2, lty=2, col=2)
    box()
    plot(uEtrains, Etrain); grid()
    abline(a=0, b=-2, lty=2, col=2)
    abline(a=0, b=2, lty=2, col=2)
    box()
    par(mfrow = c(1,1))
  }
  
  set.seed(123) # for bootstrap reprod
  resC0 = ErrViewLib::plotLZMS(
    uEtrain, Etrain/uEtrain,
    aux = 1:length(Etrain),
    logX = TRUE,
    ylim = c(0,2),
    nBin = nBin,
    plot = plot,
    popMin = 20,
    score = TRUE,
    xlab = "uE",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0('Training set, fv = ',signif(resC0$fVal,2)))

  set.seed(123)
  resC1 = ErrViewLib::plotLZMS(
    uEtrains, Etrain/uEtrains,
    aux = 1:length(Etrain),
    logX = TRUE,
    ylim = c(0,2),
    nBin = nBin,
    plot = plot,
    popMin = 20,
    score = TRUE,
    xlab = "uE scaled",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0('Scaled training set, fv = ',signif(resC1$fVal,2)))

  set.seed(123)
  resX10 = ErrViewLib::plotLZMS(
    X1train, Etrain/uEtrain,
    aux = uEtrain,
    logX = FALSE,
    ylim = c(0,2),
    nBin = nBin,
    popMin = 20,
    plot = plot,
    score = TRUE,
    xlab = "X1",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0('Training set, fv = ',signif(resX10$fVal,2)))
  
  set.seed(123)
  resX11 =ErrViewLib::plotLZMS(
    X1train, Etrain/uEtrains,
    aux = uEtrain,
    logX = FALSE,
    ylim = c(0,2),
    nBin = nBin,
    popMin = 20,
    plot = plot,
    score = TRUE,
    xlab = "X1",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0('Scaled training set, fv = ',signif(resX11$fVal,2)))

  set.seed(123)
  resX20 = ErrViewLib::plotLZMS(
    X2train, Etrain/uEtrain,
    aux = uEtrain,
    logX = FALSE,
    ylim = c(0,2),
    nBin = nBin,
    popMin = 20,
    plot = plot,
    score = TRUE,
    xlab = "X2",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0('Training set, fv = ',signif(resX20$fVal,2)))
  
  set.seed(123)
  resX21 =ErrViewLib::plotLZMS(
    X2train, Etrain/uEtrains,
    aux = uEtrain,
    logX = FALSE,
    ylim = c(0,2),
    nBin = nBin,
    popMin = 20,
    plot = plot,
    score = TRUE,
    xlab = "X2",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0('Scaled training set, fv = ',signif(resX21$fVal,2)))

  return(
    list(
      fvC0    = signif(resC0$fVal,2),
      lofvC0  = signif(resC0$lofVal,2),
      upfvC0  = signif(resC0$upfVal,2),
      fvC1    = signif(resC1$fVal,2),
      lofvC1  = signif(resC1$lofVal,2),
      upfvC1  = signif(resC1$upfVal,2),
      fvX10   = signif(resX10$fVal,2),
      lofvX10 = signif(resX10$lofVal,2),
      upfvX10 = signif(resX10$upfVal,2),
      fvX11   = signif(resX11$fVal,2),
      lofvX11 = signif(resX11$lofVal,2),
      upfvX11 = signif(resX11$upfVal,2),
      fvX20   = signif(resX20$fVal,2),
      lofvX20 = signif(resX20$lofVal,2),
      upfvX20 = signif(resX20$upfVal,2),
      fvX21   = signif(resX21$fVal,2),
      lofvX21 = signif(resX21$lofVal,2),
      upfvX21 = signif(resX21$upfVal,2)
    )
  )
}
fvScores <- function(E, u, X1, X2, nBin, method = 'stud', plot = FALSE) {
  
  Z = E / u
  M = length(Z)

  set.seed(123) # for bootstrap reprod
  resU = ErrViewLib::plotLZMS(
    u, Z,
    aux = 1:M,
    logX = TRUE,
    ylim = c(0,2),
    nBin = nBin,
    plot = plot,
    score = TRUE,
    xlab = "uE",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0(
      'Training set, fv = ',signif(resU$fVal,2),
      ' [',round(resU$lofVal,2),', ',round(resU$upfVal,2),']')
    )
  
  set.seed(123)
  resX1 = ErrViewLib::plotLZMS(
    X1, Z,
    aux = 1:M,
    logX = FALSE,
    ylim = c(0,2),
    nBin = nBin,
    plot = plot,
    score = TRUE,
    xlab = "X1",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0(
      'Training set, fv = ',signif(resX1$fVal,2),
      ' [',round(resX1$lofVal,2),', ',round(resX1$upfVal,2),']')
    )
  
  set.seed(123)
  resX2 = ErrViewLib::plotLZMS(
    X2, Z,
    aux = 1:M,
    logX = FALSE,
    ylim = c(0,2),
    nBin = nBin,
    plot = plot,
    score = TRUE,
    xlab = "X2",
    method = method,
    title = ''
  )
  if(plot)
    title(main = paste0(
      'Training set, fv = ',signif(resX2$fVal,2),
      ' [',round(resX2$lofVal,2),', ',round(resX2$upfVal,2),']')
    )
  
  return(
    list(
      fvu    = signif(resU$fVal,2),
      lofvu  = signif(resU$lofVal,2),
      upfvu  = signif(resU$upfVal,2),
      fvX1   = signif(resX1$fVal,2),
      lofvX1 = signif(resX1$lofVal,2),
      upfvX1 = signif(resX1$upfVal,2),
      fvX2   = signif(resX2$fVal,2),
      lofvX2 = signif(resX2$lofVal,2),
      upfvX2 = signif(resX2$upfVal,2)
    )
  )
}
calibrate <- function(
  Etrain, uEtrain, X1train, X2train,
  nBinOpt, intrvOpt, nBinStat, intrvStat,
  scaleMod, scoreComp, scoreFun, fv_method) {

  res = optScaleFun(
    Etrain, uEtrain, X1train, X2train,
    scaleMod, scoreComp, scoreFun,
    nBinOpt, intrvOpt, nBinStat, intrvStat
  )

  start = scoreFun(res$s0, Etrain, uEtrain, X1train, X2train, 
                   nBinOpt, intrvOpt, nBinStat, intrvStat,
                   scoreComp, scaleMod)
  end0  = scoreFun(res$s , Etrain, uEtrain, X1train, X2train, 
                   nBinOpt, intrvOpt, nBinStat, intrvStat,
                   scoreComp, scaleMod)
  print(c(start, end0, res$end))

  # Rescale train uncertainties using bins
  uEtrains = fScale(res$s, uEtrain, scaleMod, nBinOpt, intrvOpt)

  # Final scores
  
  gs   = globScore(Etrain, uEtrains,  X1train, X2train, nBinStat, intrvStat)
  fs   = fvScores(Etrain, uEtrains, X1train, X2train, nBinStat, method = 'bootstrap')
  scores = formatScores(gs, fs)
  
  return(
    list(
      scores = scores,
      s = res$s
      )
  )

}
showTabRes <- function(tabRes, Etrain, uEtrain, nBin, intrv, io0, io1, io2) {
  gs0 = globScore(Etrain, uEtrain, nBin, intrv, io0, io1, io2)
  df = data.frame(
    model = '0000',
    NLL   = signif(gs0$nll,3),
    Scal  = signif(gs0$Scal,2),
    Scon  = signif(gs0$Scon,2),
    SX1   = signif(gs0$SX1,2),
    SX2   = signif(gs0$SX2,2),
    Stot  = signif(gs0$Stot,2),
    fvuE  = paste0('[',round(tabRes[[1]]$lofvC0,2),', ',
                   round(tabRes[[1]]$upfvC0,2),']'),
    fvX1  = paste0('[',round(tabRes[[1]]$lofvX10,2),', ',
                   round(tabRes[[1]]$upfvX10,2),']'),
    fvX2  = paste0('[',round(tabRes[[1]]$lofvX20,2),', ',
                   round(tabRes[[1]]$upfvX20,2),']')
    
  )
  for(comp in names(tabRes)) {
    df1 = data.frame(
      model = comp,
      NLL   = signif(tabRes[[comp]]$nll,3),
      Scal  = signif(tabRes[[comp]]$Scal,2),
      Scon  = signif(tabRes[[comp]]$Scon,2),
      SX1   = signif(tabRes[[comp]]$SX1,2),
      SX2   = signif(tabRes[[comp]]$SX2,2),
      Stot  = signif(tabRes[[comp]]$Stot,2),
      fvuE  = paste0('[',round(tabRes[[comp]]$lofvC1,2),', ',
                     round(tabRes[[comp]]$upfvC1,2),']'),
      fvX1  = paste0('[',round(tabRes[[comp]]$lofvX11,2),', ',
                     round(tabRes[[comp]]$upfvX11,2),']'),
      fvX2  = paste0('[',round(tabRes[[comp]]$lofvX21,2),', ',
                     round(tabRes[[comp]]$upfvX21,2),']')
      
    )
    df = rbind(df,df1)
  }
  df
}
formatScores <- function(gs,fs) {
  return(
    list(
      NLL  = round(gs$nll,2),
      Scal = round(gs$Scal,2),
      Scon = round(gs$Scon,2),
      SX1  = round(gs$SX1,2),
      SX2  = round(gs$SX2,2),
      Stot = round(gs$Stot,2),
      ENCE = round(gs$ence,2),
      fvu  = paste0('[',round(fs$lofvu,2), ', ',round(fs$upfvu,2),']'),
      fvu  = paste0('[',round(fs$lofvX1,2),', ',round(fs$upfvX1,2),']'),
      fvu  = paste0('[',round(fs$lofvX2,2),', ',round(fs$upfvX2,2),']')
    )
  )  
}
fSel <- function(i, X, lwlims) {
  if(i == 1)
    xmin = -Inf
  else
    xmin = lwlims[i-1]

  if(i == length(lwlims) +1 )
    xmax  = Inf
  else
    xmax = lwlims[i]

  sel = X >= xmin & X < xmax

  if(sum(sel) == 0)
    print(c(i, xmin, xmax))

  return(sel)
}
BVSfactors <- function(X, Z, lwlims, nb) {
  s = c()
  for (i in 1:nb) {
    sel  = fSel(i, X, lwlims)
    s[i] = sqrt(mean(Z[sel]^2))
  }
  return(s)
}

# Get training and test sets ####
D = read.table(
  '../Data/BUS2022/qm9_qm9_E_uncalibrated_validation.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

S  = D[, "formula"]
X1 = CHNOSZ::mass(S)
X2 = fractHetero(S)
R  = D[, "E_"]
C  = D[, "prediction"]
uC = D[, "uncertainty_total"] ^ 0.5
E  = R - C
M = length(E)
uE = uC

D = read.table(
  '../Data/BUS2022/qm9_qm9_E_uncalibrated_test.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

St  = D[, "formula"]
X1t = CHNOSZ::mass(St)
X2t = fractHetero(St)
Rt  = D[, "E_"]
Ct  = D[, "prediction"]
uCt = D[, "uncertainty_total"] ^ 0.5
Et  = Rt - Ct
Mt = length(Et)
uEt = uCt

D = read.table(
  '../Data/BUS2022/qm9_qm9_E_calibrated_isotonic_test.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

Si  = D[, "formula"]
X1i = CHNOSZ::mass(Si)
X2i = fractHetero(Si)
Ri  = D[, "E_"]
Ci  = D[, "prediction"]
uCi = D[, "uncertainty_total"] ^ 0.5
Ei  = Ri - Ci
uEi = uCi

# Scoring bins
nBinStat = 100
intrvStat = ErrViewLib::genIntervals(M, nBinStat)

fv_method = 'bootstrap'
# fv_method = 'stud'

scaleMod= 'bin-wise'

nMC = 1000

# uE binning ####

io = order(uE)
uEtrain = uE[io]
Etrain  = E[io]
X1train = X1[io]
X2train = X2[io]

io = order(uEt)
uEtest = uEt[io]
Etest  = Et[io]
X1test = X1t[io]
X2test = X2t[io]

io    = order(uEi)
uEiso = uEi[io]
Eiso  = Ei[io]
X1iso = X1i[io]
X2iso = X2i[io]

## Reference values ####
tabResAux = list()
gs = globScore(Etrain, uEtrain, X1train, X2train, nBinStat, intrvStat)
fs = fvScores(Etrain, uEtrain, X1train, X2train, nBinStat, method = fv_method)
tabResAux[['train_0000']] = formatScores(gs, fs)

gs = globScore(Etest, uEtest, X1test, X2test, nBinStat, intrvStat)
fs = fvScores(Etest, uEtest, X1test, X2test, nBinStat, method = fv_method)
tabResAux[['test_0000']] = formatScores(gs, fs)

gs = globScore(Eiso, uEiso, X1iso, X2iso, nBinStat, intrvStat)
fs = fvScores(Eiso, uEiso, X1iso, X2iso, nBinStat, method = fv_method)
tabResAux[['iso_0000']] = formatScores(gs, fs)

tabSim = list()
for (i in 1:nMC) {
  Esim = rnorm(uEtrain, 0, uEtrain)
  gs = globScore(Esim, uEtrain, X1train, X2train, nBinStat, intrvStat)
  tabSim[[i]]  = gs
}
tab = as.data.frame(do.call(rbind, tabSim))
gs = as.list(apply(tab,2,function(x) mean(unlist(x))))
fs = fvScores(Esim, uEtrain, X1train, X2train, nBinStat, method = 'stud')
tabResAux[['train_simul_0000']] = formatScores(gs, fs)

tabSim = list()
for (i in 1:nMC) {
  Esim = rnorm(uEtest,0,uEtest)
  gs = globScore(Esim, uEtest, X1test, X2test, nBinStat, intrvStat)
  tabSim[[i]]  = gs
}
tab = as.data.frame(do.call(rbind, tabSim))
gs = as.list(apply(tab,2,function(x) mean(unlist(x))))
fs = fvScores(Esim, uEtest, X1test, X2test, nBinStat, method = 'stud')
tabResAux[['test_simul_0000']] = formatScores(gs, fs)

tab = as.matrix(as.data.frame(do.call(rbind, tabResAux)))
knitr::kable(tab)

models  = c("2000","1111","1100","1011")

for(nBinOpt in c(20,40,80)) {
  intrvOpt = ErrViewLib::genIntervals(uEtrain, nBinOpt)
  # nBinOpt  = intrvOpt$nbr # to use when equiPop = FALSE

  tabRes = tabResAux 
  tabBVS = list()
  for(scoreComp in models) {
    
    # Optimize scaling factors
    res = calibrate(
      Etrain, uEtrain, X1train, X2train,
      nBinOpt, intrvOpt, nBinStat, intrvStat,
      scaleMod, scoreComp, scoreFun, fv_method)
    tabRes[[paste0('train_',scoreComp)]] = res$scores
    tabBVS[[scoreComp]] = res$s
    
    tab = as.matrix(as.data.frame(do.call(rbind, tabRes)))
    knitr::kable(tab)
  } # Separate loops to order results table rows
  
  for(scoreComp in models) {
    # Scale test set using training set intervals limits
    s = tabBVS[[scoreComp]]
    lwlims  = uEtrain[intrvOpt$lwindx][-1]
    uEtests = scaleByLimits(uEtest, uEtest, s, lwlims)
    
    # Score scaled test set
    gs = globScore(Etest, uEtests, X1test, X2test, nBinStat, intrvStat)
    fs = fvScores( Etest, uEtests, X1test, X2test, nBinStat, method = fv_method)
    tabRes[[paste0('test_',scoreComp)]] = formatScores(gs, fs)
  
    tab = as.matrix(as.data.frame(do.call(rbind, tabRes)))
    knitr::kable(tab)
  }

  sink(file =  file.path(tabDir,paste0('tabRes_',nBinOpt,'.tex')))
  tab = as.matrix(as.data.frame(do.call(rbind, tabRes)))
  print(knitr::kable(tab, 'latex'))
  sink()
}

save(tabRes,tabBVS,file=paste0(tabDir,'/tab_Res_BVS.Rda'))


# X1 binning ####
io = order(X1,uE)
uEtrain = uE[io]
Etrain  = E[io]
X1train = X1[io]
X2train = X2[io]

io = order(X1t,uEt)
uEtest = uEt[io]
Etest  = Et[io]
X1test = X1t[io]
X2test = X2t[io]

models2 =  c("2000","1111","1100","1011")

for(nBinOpt in c(20,40,80)) {
  scaleMod = 'bin-wise'
  intrvOpt  = ErrViewLib::genIntervals(uEtrain, nBinOpt)
  # nBinOpt = intrvOpt$nbr
  
  tabRes2 = tabResAux 
  tabBVS2 = list()
  for(scoreComp in models2) {
    
    # Optimize scaling factors
    res = calibrate(
      Etrain, uEtrain, X1train, X2train,
      nBinOpt, intrvOpt, nBinStat, intrvStat,
      scaleMod, scoreComp, scoreFun, fv_method)
    tabRes2[[paste0('train_',scoreComp)]] = res$scores
    tabBVS2[[scoreComp]] = res$s
    
    tab = as.matrix(as.data.frame(do.call(rbind, tabRes2)))
    knitr::kable(tab)
  } # Separate loops to order results table rows
  
  
  # Rescale test uncertainties using X1train and uEtrain limits
  io = order(X1,uE)
  X1train = X1[io]
  uEtrain = uE[io]
  for(scoreComp in models2) {
    
    s = tabBVS2[[scoreComp]]
    
    lwlims1 = X1train[intrvOpt$lwindx][-1]
    lwlims2 = uEtrain[intrvOpt$lwindx][-1]
    uEtests = uEtest
    for(i in 1:nBinOpt) {
      if(i == 1)
        xmin = -Inf
      else
        xmin = lwlims1[i-1]
      
      if(i == nBinOpt)
        xmax  = Inf
      else
        xmax = lwlims1[i]
      
      if(xmin < xmax) {
        sel = X1test >= xmin & X1test < xmax
        
      } else { # degeneracy due to stratification: use secondary variable
        sel = X1test == xmin
        
        if (sum(lwlims1 == xmin) > 2) {
          
          indxs = which(lwlims1 == xmin)
          
          if(i == indxs[1])
            xmin2 = -Inf
          else
            xmin2 = lwlims2[i-1]
          
          if(i == indxs[length(indxs)])
            xmax2  = Inf
          else
            xmax2 = lwlims2[i]
          
          sel = sel & uEtest >= xmin2 & uEtest < xmax2
          
        }
      }
      if(sum(sel) > 0)
        uEtests[sel] = uEtest[sel] * s[i]
      
    }
    
    gs = globScore(Etest, uEtests, X1test, X2test, nBinStat, intrvStat)
    fs = fvScores( Etest, uEtests, X1test, X2test, nBinStat, method = fv_method)
    tabRes2[[paste0('test_',scoreComp)]] = formatScores(gs, fs)
    
    tab = as.matrix(as.data.frame(do.call(rbind, tabRes2)))
    knitr::kable(tab)
  }
  
  sink(file =  file.path(tabDir,paste0('tabRes2_',nBinOpt,'.tex')))
  tab = as.matrix(as.data.frame(do.call(rbind, tabRes2)))
  print(knitr::kable(tab, 'latex'))
  sink()
  
}
save(tabRes2,tabBVS2,file=paste0(tabDir,'/tab_Res2_BVS2.Rda'))

# Fig 1 ####
png(
  file = file.path(figDir, 'Fig_01a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E ,uE ,
  score = FALSE,
  probref = FALSE,
  showLegend = FALSE,
  unit = '[eV]',
  label = 1,
  title = 'Confidence curve',
  gPars = gPars
)
ErrViewLib::plotCC(
  Et, uEt,
  score = FALSE,
  probref = FALSE,
  add = TRUE,
  col =2,
  gPars = gPars
)
legend('topright', bty ='n',
       legend = c('Train','Test'),
       col = gPars$cols[c(6,2)],
       lty =1, lwd = gPars$lwd
)
dev.off()

png(
  file = file.path(figDir, 'Fig_01b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLZMS(
  uEtrain, Etrain/uEtrain,
  logX = TRUE,
  xlab = 'Uncertainty, u [eV]',
  nBin = nBinStat,
  score = FALSE,
  label = 2,
  title = 'LZMS analysis',
  gPars = gPars
)
dev.off()

png(
  file = file.path(figDir, 'Fig_01c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLZMS(
  X1train, Etrain/uEtrain,
  logX = FALSE,
  xlab = 'X1, [Da]',
  nBin = nBinStat,
  score = FALSE,
  label = 3,
  title = 'LZMS analysis',
  gPars = gPars
)
dev.off()

png(
  file = file.path(figDir, 'Fig_01d.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLZMS(
  X2train, Etrain/uEtrain,
  logX = FALSE,
  xlab = 'X2',
  nBin = nBinStat,
  score = FALSE,
  label = 4,
  title = 'LZMS analysis',
  gPars = gPars
)
dev.off()

# Optimal binning ####

io = order(uE)
uEtrain = uE[io]
Etrain  = E[io]
X1train = X1[io]
X2train = X2[io]

io = order(uEt)
uEtest = uEt[io]
Etest  = Et[io]
X1test = X1t[io]
X2test = X2t[io]

nBins = c(2, 3, 4, 5, seq(6, 80, by = 2))
co = matrix(NA, ncol = 19, nrow = length(nBins))
nb = c()
# fv_method = 'stud'
nMC    = 1000

for (j in seq_along(nBins)) {
  nBinOpt  = nBins[j]; print(nBinOpt)
  intrvOpt = ErrViewLib::genIntervals(uEtrain, nBinOpt)
  lwlims   = uEtrain[intrvOpt$lwindx][-1]
  nb[j]    = intrvOpt$nbr
  s        = BVSfactors(uEtrain, Etrain / uEtrain, lwlims, nb[j])
  uEtests  = scaleByLimits(uEtest, uEtest, s, lwlims)
  gs       = globScore(Etest, uEtests, X1test, X2test, nBinStat, intrvStat)
  fs       = fvScores( Etest, uEtests, X1test, X2test, nBinStat, method = fv_method)
  co[j, 1:10] = unlist(gs)
  co[j,11:19] = unlist(fs)
}
save(co, file = paste0(tabDir,'/tab_co.Rda'))

## Simulated datasets ####
nBins = c(2, 3, 4, 5, seq(6, 80, by = 2))
cosim = matrix(NA, ncol = 10, nrow = length(nBins))
nb = c()
nMC = 10000
for (j in seq_along(nBins)) {
  nBinOpt  = nBins[j]; print(nBinOpt)
  intrvOpt = ErrViewLib::genIntervals(uEtrain, nBinOpt)
  lwlims   = uEtrain[intrvOpt$lwindx][-1]
  nb[j]    = intrvOpt$nbr
  s        = BVSfactors(uEtrain, Etrain / uEtrain, lwlims, nb[j])
  uEtests  = scaleByLimits(uEtest, uEtest, s, lwlims)

  tabSim = list()
  for (i in 1:nMC) {
    Esim = rnorm(uEtests,0,uEtests)
    gs   = globScore(Esim, uEtests, X1test, X2test, nBinStat, intrvStat)
    tabSim[[i]]  = gs
  }
  tab = as.data.frame(do.call(rbind, tabSim))
  gs  = apply(tab,2,function(x) mean(unlist(x)))
  cosim[j, ] = gs
}
save(cosim, file = paste0(tabDir,'/tab_cosim.Rda'))

# Fig 2 ####

png(
  file = file.path(figDir, paste0('Fig_02.png')),
  width  = 2 * gPars$reso,
  height = 2 * gPars$reso
)
par(
  mfrow = c(2,2),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)

## Fig2a ####
matplot( #NLL
  nb,
  co[, 8],
  type = 'p',
  pch = 19,
  xlab = 'Nb. of bins',
  ylab = 'Negative log likelihood',
  ylim = c(-3.08,-3.01),
  main = 'NLL',
  col = gPars$cols[1]
)
grid()
abline(h = tabResIso$iso$nll,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[1])
matlines(nb, cosim[,8], lty = 3, lwd =  3* gPars$lwd, col = gPars$cols[1])
legend(
  'topright',
  bty    = 'n',
  legend = c('Score', 'Isotonic', 'Simul'),
  pch    = c(19, NA, NA),
  lty    = c(0, 2, 3),
  lwd    = c(0, 2, 3) * gPars$lwd,
  col    = gPars$cols[c(1, 1, 1)]
)
box()
label = 1
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.25
)

## Fig2b ####
matplot(
  nb,
  cbind(co[, 6], 1e3*co[,7]),
  type = 'p',
  pch = 19,
  col = gPars$cols[1:2],
  xlab = 'Nb. of bins',
  ylab = 'Calibration errors',
  ylim = c(0, 0.15),
  yaxs = 'i',
  main = 'ENCE, UCE'
)
grid()
abline(h = tabResIso$iso$ence,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[1])
abline(h = tabResIso$iso$uce*1e3,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[2])
matlines(
  nb,
  cbind(cosim[, 6], 1e3*cosim[,7]),
  type = 'l',
  lty = 3, lwd = 3* gPars$lwd,
  col = gPars$cols[1:2]
)
legend(
  'top',
  bty    = 'n',
  ncol   = 2,
  legend = c('ENCE', 'UCE * 1000'),
  pch    = 19,
  lty    = 0,
  col    = gPars$cols[1:2]
)
box()
label = 2
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.25
)

## Fig2c ####
matplot(
  nb,
  cbind(co[, 1:4],co[,5]/2),
  type = 'p',
  pch = 19,
  col = gPars$cols[1:5],
  xlab = 'Nb. of bins',
  ylab = 'Calibration scores',
  ylim = c(0, 0.6),
  yaxs = 'i',
  main = expression(S[x])
)
grid()
abline(h = tabResIso$iso$Scal,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[1])
abline(h = tabResIso$iso$Scon,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[2])
abline(h = tabResIso$iso$SX1,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[3])
abline(h = tabResIso$iso$SX2,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[4])
abline(h = tabResIso$iso$Stot/2,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[5])
matlines(
  nb,
  cbind(cosim[, 1:4],cosim[,5]/2),
  type = 'l',
  lty = 3, lwd = 3* gPars$lwd,
  col = gPars$cols[1:5]
)
legend(
  'topright',
  bty = 'n',
  ncol = 3,
  legend = c(
    expression(S[cal]),
    expression(S[u]),
    expression(S[X[1]]),
    expression(S[X[2]]),
    expression(S[tot]/2)
  ),
  pch = 19,
  lty = 0,
  col = gPars$cols[1:5]
)
box()
label = 3
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.25
)

## Fig2d ####
plot(
  nb,
  co[, 11],
  type = 'p',
  pch = 19,
  col = gPars$cols[2],
  xlab = 'Nb. of bins',
  ylab = 'Percentage of valid intervals',
  ylim = c(0.5, 1.066),
  yaxs = 'i',
  main = expression(f[list(v, x)])
)
grid()
segments(nb, co[, 12], nb, co[, 13], col = gPars$cols[2])
abline(h = tabResIso$iso$fvC1,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[2])
points(nb,
       co[, 14],
       type = 'p',
       pch = 19,
       col = gPars$cols[3])
segments(nb, co[, 15], nb, co[, 16], col = gPars$cols[3])
abline(h = tabResIso$iso$fvX11,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[3])
points(nb,
       co[, 17],
       type = 'p',
       pch = 19,
       col = gPars$cols[4])
segments(nb, co[, 18], nb, co[, 19], col = gPars$cols[4])
abline(h = tabResIso$iso$fvX21,
       lty = 2, lwd = 2* gPars$lwd,
       col = gPars$cols[4])
abline(h = 0.95,
       lty = 3,lwd = 3* gPars$lwd,
       col = gPars$cols[4])
polygon(c(-10,100,100,-10),
        c(0.9,0.9,1.0,1.0),
        col = gPars$cols_tr[4], border = NA)
legend(
  'top',
  bty = 'n',
  ncol = 3,
  legend = c(expression(f[list(v, u)]),
             expression(f[list(v, X[1])]),
             expression(f[list(v, X[2])])),
  pch = 19,
  lty = 0,
  col = gPars$cols[2:4]
)
box()
label = 4
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.25
)
dev.off()


# 2D binning ####
tabRes2D = list()

nBinOpt = 20
intrvOptU  = ErrViewLib::genIntervals(uEtrain, nBinOpt)
nBinOptu = intrvOptU$nbr
intrvOptX  = ErrViewLib::genIntervals(X1train, nBinOpt)
nBinOptX = intrvOptX$nbr

# Training
savg    = sqrt(mean((Etrain/uEtrain)^2))
lwlimsU = sort(uEtrain)[intrvOptU$lwindx][-1]
lwlimsX = sort(X1train)[intrvOptX$lwindx][-1]
s0 = matrix(NA,nrow = nBinOptu, ncol = nBinOptX)
for(i in 1:nBinOptu) {
  seli  = fSel(i, uEtrain, lwlimsU)
  for(j in 1:nBinOptX) {
    selj  = fSel(j, X1train, lwlimsX)
    sel = seli & selj
    if(sum(sel) == 0)
      s0[i,j] = savg
    else
      s0[i,j] = sqrt(mean((Etrain[sel]/uEtrain[sel])^2))
  }
}

# Score
uEtrains = uEtrain
for(i in 1:nBinOptu) {
  seli  = fSel(i, uEtrain, lwlimsU)
  for(j in 1:nBinOptX) {
    selj  = fSel(j, X1train, lwlimsX)
    sel = seli & selj
    if(sum(sel) != 0)
      uEtrains[sel] = uEtrain[sel] * s0[i,j]
  }
}

gs   = globScore(Etrain, uEtrains, X1train, X2train, nBinStat, intrvStat)
fs   = fvScores(Etrain, uEtrains, X1train, X2train, nBinStat, method = fv_method)
tabRes2D[['train']] = formatScores(gs, fs)

# Test
uEtests = uEtest
for(i in 1:nBinOptu) {
  seli  = fSel(i, uEtest, lwlimsU)
  for(j in 1:nBinOptX) {
    selj  = fSel(j, X1test, lwlimsX)
    sel = seli & selj
    if(sum(sel) != 0)
      uEtests[sel] = uEtest[sel] * s0[i,j]
  }
}
gs = globScore(Etest, uEtests, X1test, X2test, nBinStat, intrvStat)
fs = fvScores( Etest, uEtests, X1test, X2test, nBinStat, method = fv_method)
tabRes2D[['test']] = formatScores(gs, fs)

sink(file =  file.path(tabDir,paste0('tabRes2D.tex')))
tab = as.matrix(as.data.frame(do.call(rbind, tabRes2D)))
print(knitr::kable(tab, 'latex'))
sink()

# resp = plotRes(Etest, uEtest, uEtests, X1test, X2test, 
#                nBinStat, method = 'stud', plot = TRUE)
# fvScores(Etrain, uEtrain, X1train, X2train,
#          nBinStat, method = 'bootstrap', plot=TRUE)
# fvScores(Etrain, uEtrains, X1train, X2train,
#          nBinStat, method = 'bootstrap', plot=TRUE)
# fvScores(Etest, uEtests, X1test, X2test,
#          nBinStat, method = 'bootstrap', plot=TRUE)
