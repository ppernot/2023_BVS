
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
      ence = mean(abs(sqrt(MV)-sqrt(MSE))/sqrt(MV)),
      uce  = mean(abs(MV-MSE))
    )
  )
}
nll = function(E, uE) {
  0.5*mean( E^2/uE^2 + log(uE^2) + 1.837877 )
}
scoreFun = function(
  pars,
  E, uE,
  nBinOpt, intrvOpt, nBinStat, intrvStat,
  io0, io1, io2,
  scoreComp ='1111',
  scaleMod = c('bin-wise','element-wise','linear')) {

  scaleMod = match.arg(scaleMod)

  # Apply scaling model to uE
  uE   = fScale(pars, uE, scaleMod, nBinOpt, intrvOpt)

  score = 0

  # Mean calibration error
  if (substr(scoreComp, 1, 1) == "1")
    score = score + abs(log(mean((E / uE) ^ 2)))
  else if (substr(scoreComp, 1, 1) == "2")
    score = score + nll(E, uE)

  # Mean consistency error
  if (substr(scoreComp, 2, 2) == "1")
    score = score + scoreFun0(E[io0] / uE[io0], nBinStat, intrvStat)

  # Mean adaptivity error
  ## X1
  if (substr(scoreComp, 3, 3) == "1")
    score = score + scoreFun0(E[io1] / uE[io1], nBinStat, intrvStat)
  ## X2
  if (substr(scoreComp, 4, 4) == "1")
    score = score + scoreFun0(E[io2] / uE[io2], nBinStat, intrvStat)

  return(score)
}
globScore <- function(E, uE, nBin, intrv, io0, io1, io2) {
  Z = E / uE
  cal   = abs(log(mean(Z^2)))
  cons  = scoreFun0(Z[io0], nBin, intrv)
  adap1 = scoreFun0(Z[io1], nBin, intrv)
  adap2 = scoreFun0(Z[io2], nBin, intrv)
  li    = scoreFun1(E[io0],uE[io0], nBin, intrv)
  ence  = li$ence
  uce   = li$uce
  return(
    c(cal, cons, adap1, adap2, nll(E,uE), ence, uce)
  )
}
optScaleFun <- function(
  Etrain, uEtrain,
  scaleMod, scoreComp, scoreFun,
  nBinOpt, intrvOpt, nBinStat, intrvStat,
  io0, io1, io2
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
    start = scoreFun(s0, Etrain, uEtrain,
                     nBinOpt, intrvOpt, nBinStat, intrvStat,
                     io0, io1, io2, scoreComp, scaleMod)

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
    E = Etrain, uE = uEtrain,
    nBinOpt = nBinOpt, intrvOpt = intrvOpt,
    nBinStat = nBinStat, intrvStat = intrvStat,
    io0 = io0, io1 = io1, io2 = io2,
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
      E = Etrain, uE = uEtrain,
      nBinOpt = nBinOpt, intrvOpt = intrvOpt,
      nBinStat = nBinStat, intrvStat = intrvStat,
      io0 = io0, io1 = io1, io2 = io2,
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
                    method = 'bootstrap', plot = FALSE) {

  # Control plots, but mostly to generate ZMS validation stats

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
  title(main = paste0('Scaled training set, fv = ',signif(resC1$fVal,2)))

  set.seed(123)
  resX10 = ErrViewLib::plotLZMS(
    X1train, Etrain/uEtrain,
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
  title(main = paste0('Training set, fv = ',signif(resX10$fVal,2)))
  set.seed(123)
  resX11 =ErrViewLib::plotLZMS(
    X1train, Etrain/uEtrains,
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
  title(main = paste0('Scaled training set, fv = ',signif(resX11$fVal,2)))
  set.seed(123)
  resX20 = ErrViewLib::plotLZMS(
    X2train, Etrain/uEtrain,
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
  title(main = paste0('Training set, fv = ',signif(resX20$fVal,2)))
  set.seed(123)
  resX21 =ErrViewLib::plotLZMS(
    X2train, Etrain/uEtrains,
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
calibrate <- function(
  Etrain, uEtrain, X1train, X2train,
  nBinOpt, intrvOpt, nBinStat, intrvStat,
  io0, io1, io2,
  scaleMod, scoreComp, scoreFun, tabRes, tabBVS) {

  res = optScaleFun(
    Etrain, uEtrain,
    scaleMod, scoreComp, scoreFun,
    nBinOpt, intrvOpt, nBinStat, intrvStat,
    io0, io1, io2
  )

  start = scoreFun(res$s0, Etrain, uEtrain,
                   nBinOpt, intrvOpt, nBinStat, intrvStat,
                   io0, io1, io2, scoreComp, scaleMod)
  end0  = scoreFun(res$s , Etrain, uEtrain,
                   nBinOpt, intrvOpt, nBinStat, intrvStat,
                   io0, io1, io2, scoreComp, scaleMod)
  print(c(start, end0, res$end))
  plot(res$s0/res$s); abline(h=1,lty=2,col=2)
  tabBVS[[scoreComp]] = res$s

  # Rescale train uncertainties using bins
  uEtrains = fScale(res$s, uEtrain, scaleMod, nBinOpt, intrvOpt)
  plot(uEtrain, uEtrains, log = 'xy', col = 4)
  abline(a=0,b=1, lty = 2)
  grid(equilogs = FALSE)

  # Final sub-scores
  abs(log(mean((Etrain/uEtrains)^2)))
  scoreFun0(Etrain/uEtrains, nBinStat, intrvStat)
  scoreFun0(Etrain[io1]/uEtrains[io1], nBinStat, intrvStat)
  scoreFun0(Etrain[io2]/uEtrains[io2], nBinStat, intrvStat)

  resp = plotRes(Etrain, uEtrain, uEtrains, X1train, X2train, nBinStat)
  gs   = globScore(Etrain, uEtrains, nBinStat, intrvStat, io0, io1, io2)
  resp$cal   = gs[1]
  resp$cons  = gs[2]
  resp$adap1 = gs[3]
  resp$adap2 = gs[4]
  resp$glob  = sum(gs[1:4])
  resp$nll   = gs[5]
  tabRes[[scoreComp]] = resp

  return(
    list(
      tabRes = tabRes,
      tabBVS = tabBVS
      )
  )

}
showTabRes <- function(tabRes, Etrain, uEtrain, nBin, intrv, io0, io1, io2) {
  gs0 = globScore(Etrain, uEtrain, nBin, intrv, io0, io1, io2)
  df = data.frame(
    model = '0000',
    nll   = signif(gs0[5],3),
    cal   = signif(gs0[1],2),
    cons  = signif(gs0[2],2),
    adap1 = signif(gs0[3],2),
    adap2 = signif(gs0[4],2),
    glob  = signif(sum(gs0[1:4]),2),
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

      nll   = signif(tabRes[[comp]]$nll,3),
      cal   = signif(tabRes[[comp]]$cal,2),
      cons  = signif(tabRes[[comp]]$cons,2),
      adap1 = signif(tabRes[[comp]]$adap1,2),
      adap2 = signif(tabRes[[comp]]$adap2,2),
      glob  = signif(tabRes[[comp]]$glob,2),
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

# uE binning ####

io = order(uE)
uEtrain = uE[io]
Etrain  = E[io]
X1train = X1[io]
X2train = X2[io]
io0 = order(uEtrain)
io1 = order(X1train,uEtrain)
io2 = order(X2train,uEtrain)

io = order(uEt)
uEtest = uEt[io]
Etest  = Et[io]
X1test = X1t[io]
X2test = X2t[io]
io0t = order(uEtest)
io1t = order(X1test,uEtest)
io2t = order(X2test,uEtest)

nBinStat = 100
intrvStat = ErrViewLib::genIntervals(M, nBinStat)
models = c("2000","1111","1100","1011")

for(nBinOpt in c(20,40,80)) {
  scaleMod = 'bin-wise'
  intrvOpt  = ErrViewLib::genIntervals(uEtrain, nBinOpt,
                                       equiPop = TRUE,
                                       logBin = TRUE)
  nBinOpt = intrvOpt$nbr

  tabRes = tabBVS = list()
  for(scoreComp in models) {
    res = calibrate(
      Etrain, uEtrain, X1train, X2train,
      nBinOpt, intrvOpt, nBinStat, intrvStat,
      io0, io1, io2,
      scaleMod, scoreComp, scoreFun,
      tabRes, tabBVS)
    tabRes = res$tabRes
    tabBVS = res$tabBVS
    print(showTabRes(tabRes, Etrain, uEtrain,
                     nBinStat, intrvStat,
                     io0, io1, io2))
  }

  for(scoreComp in models) {
    s = tabBVS[[scoreComp]]
    lwlims  = uEtrain[intrvOpt$lwindx][-1]
    uEtests = scaleByLimits(uEtest, uEtest, s, lwlims)
    resp = plotRes(Etest, uEtest, uEtests, X1test, X2test, nBinStat)
    gs   = globScore(Etest, uEtests, nBinStat, intrvStat, io0t, io1t, io2t)
    resp$cal   = gs[1]
    resp$cons  = gs[2]
    resp$adap1 = gs[3]
    resp$adap2 = gs[4]
    resp$glob  = sum(gs[1:4])
    resp$nll   = gs[5]
    tabRes[[paste0('test_',scoreComp)]] = resp
  }

  sink(file =  file.path(tabDir,paste0('tabRes_',nBinOpt,'.tex')))
  tab = showTabRes(tabRes, Etrain, uEtrain, nBinStat, intrvStat, io0, io1, io2)
  print(knitr::kable(tab, 'latex'))
  sink()
}

# X1 binning ####
io = order(X1,uE)
uEtrain = uE[io]
Etrain  = E[io]
X1train = X1[io]
X2train = X2[io]
io0 = order(uEtrain)
io1 = order(X1train)
io2 = order(X2train,uEtrain)

io = order(X1t,uEt)
uEtest = uEt[io]
Etest  = Et[io]
X1test = X1t[io]
X2test = X2t[io]
io0t = order(uEtest)
io1t = order(X1test)
io2t = order(X2test,uEtest)

for(nBinOpt in c(20,40,80)) {
  scaleMod = 'bin-wise'
  intrvOpt  = ErrViewLib::genIntervals(uEtrain, nBinOpt,
                                       equiPop = TRUE,
                                       logBin = FALSE)
  nBinOpt = intrvOpt$nbr

  models2 =  c("2000","1111","1100","1011")

  tabRes2 = tabBVS2 = list()
  for(scoreComp in models2) {
    res = calibrate(
      Etrain, uEtrain, X1train, X2train,
      nBinOpt, intrvOpt, nBinStat, intrvStat,
      io0, io1, io2,
      scaleMod, scoreComp, scoreFun,
      tabRes2, tabBVS2)
    tabRes2 = res$tabRes
    tabBVS2 = res$tabBVS
    print(showTabRes(tabRes2, Etrain, uEtrain,
                     nBinStat, intrvStat,
                     io0, io1, io2))
  }

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

    resp = plotRes(Etest, uEtest, uEtests, X1test, X2test, nBinStat)
    gs   = globScore(Etest, uEtests, nBinStat, intrvStat, io0t, io1t, io2t)
    resp$cal   = gs[1]
    resp$cons  = gs[2]
    resp$adap1 = gs[3]
    resp$adap2 = gs[4]
    resp$glob  = sum(gs[1:4])
    resp$nll   = gs[5]
    tabRes2[[paste0('test_',scoreComp)]] = resp
  }

  sink(file =  file.path(tabDir,paste0('tabRes2_',nBinOpt,'.tex')))
  tab = showTabRes(tabRes2, Etrain, uEtrain, nBinStat, intrvStat, io0, io1, io2)
  print(knitr::kable(tab, 'latex'))
  sink()
}


# Isotonic reg. ####

D = read.table(
  '../Data/BUS2022/Busk2022_QM9_E.csv', #'qm9_U0_test_Orig.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

Si  = D[, "formula"]
X1i = CHNOSZ::mass(Si)
X2i = fractHetero(Si)
Ri  = D[, "E"]
Ci  = D[, "prediction"]
uCi = D[, "uncertainty_total"] ^ 0.5
Ei = Ri - Ci
uEi = uCi

io = order(uEi)
uEiso = uEi[io]
Eiso  = Ei[io]
X1iso = X1i[io]
X2iso = X2i[io]
io0i = order(uEiso)
io1i = order(X1iso,uEiso)
io2i = order(X2iso,uEiso)

resp = plotRes(Eiso, uEiso, uEiso, X1iso, X2iso, nBinStat)
gs   = globScore(Eiso, uEiso, nBinStat, intrvStat, io0i, io1i, io2i)
resp$cal   = gs[1]
resp$cons  = gs[2]
resp$adap1 = gs[3]
resp$adap2 = gs[4]
resp$glob  = sum(gs[1:4])
resp$nll   = gs[5]
resp$ence  = gs[6]
resp$uce   = gs[7]

tabResIso = list()
tabResIso[['2000']] = tabRes[['2000']]
tabResIso[['iso']]  = resp

sink(file =  file.path(tabDir,paste0('tabResIso.tex')))
tab = showTabRes(tabResIso, Etrain, uEtrain, nBinStat, intrvStat, io0, io1, io2)
print(knitr::kable(tab, 'latex'))
sink()


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

# Optimal binning ####
nBins = c(2, 3, 4, 5, seq(6, 80, by = 2))
co = matrix(NA, ncol = 16, nrow = length(nBins))
nb = c()
uEtrain = sort(uEtrain)
for (j in seq_along(nBins)) {
  nBinOpt  = nBins[j]
  intrvOpt = ErrViewLib::genIntervals(uEtrain, nBinOpt)
  lwlims   = uEtrain[intrvOpt$lwindx][-1]
  nb[j]    = intrvOpt$nbr
  s        = BVSfactors(uEtrain, Etrain / uEtrain, lwlims, nb[j])
  uEtests  = scaleByLimits(uEtest, uEtest, s, lwlims)
  gs       = globScore(Etest, uEtests,
                       nBinStat, intrvStat,
                       io0t, io1t, io2t)
  co[j, 1:7] = gs
  set.seed(123)
  res      = ErrViewLib::plotLZMS(
    uEtests,
    Etest / uEtests,
    nBin = nBinStat,
    score = TRUE,
    method = 'bootstrap',
    plot = FALSE
  )
  co[j, 8:10]  = c(res$fVal, res$lofVal, res$upfVal)
  set.seed(123)
  res      = ErrViewLib::plotLZMS(
    X1test,
    Etest / uEtests,
    nBin = nBinStat,
    score = TRUE,
    method = 'bootstrap',
    plot = FALSE
  )
  co[j, 11:13] = c(res$fVal, res$lofVal, res$upfVal)
  set.seed(123)
  res      = ErrViewLib::plotLZMS(
    X2test,
    Etest / uEtests,
    nBin = nBinStat,
    score = TRUE,
    method = 'bootstrap',
    plot = FALSE
  )
  co[j, 14:16] = c(res$fVal, res$lofVal, res$upfVal)
}

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
plot(
  nb,
  co[, 5],
  type = 'p',
  pch = 19,
  xlab = 'Nb. of bins',
  ylab = 'Negative log likelihood',
  ylim = c(-3.08,-3.01),
  main = 'NLL',
  col = gPars$cols[5]
)
grid()
abline(h = tabResIso$nll,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[5])
box()
label = 1
mtext(
  text = paste0('(', letters[label], ')'),
  side = 3,
  adj = 1,
  cex = gPars$cex,
  line = 0.25
)

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
abline(h = tabResIso$ence,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[1])
abline(h = tabResIso$uce*1e3,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[2])
legend(
  'top',
  bty = 'n',
  ncol = 2,
  legend = c('ENCE', 'UCE x 1000'),
  pch = 19,
  lty = 0,
  col = gPars$cols[1:2]
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

matplot(
  nb,
  co[, 1:4],
  type = 'p',
  pch = 19,
  col = gPars$cols[1:4],
  xlab = 'Nb. of bins',
  ylab = 'Calibration scores',
  ylim = c(0, 0.45),
  yaxs = 'i',
  main = expression(S[x])
)
grid()
abline(h = tabResIso$cal,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[1])
abline(h = tabResIso$cons,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[2])
abline(h = tabResIso$adap1,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[3])
abline(h = tabResIso$adap2,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[4])
legend(
  'top',
  bty = 'n',
  ncol = 4,
  legend = c(
    expression(S[cal]),
    expression(S[u]),
    expression(S[X[1]]),
    expression(S[X[2]])
  ),
  pch = 19,
  lty = 0,
  col = gPars$cols[1:4]
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

plot(
  nb,
  co[, 8],
  type = 'p',
  pch = 19,
  col = gPars$cols[2],
  xlab = 'Nb. of bins',
  ylab = 'Percentage of valid intervals',
  ylim = c(0.5, 1),
  yaxs = 'i',
  main = expression(f[list(v, x)])
)
grid()
segments(nb, co[, 9], nb, co[, 10], col = gPars$cols[2])
abline(h = tabResIso$fvC1,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[2])
polygon(c(-10,100,100,-10),
        c(tabResIso$lofvC1,tabResIso$lofvC1,
          tabResIso$upfvC1,tabResIso$upfvC1),
        col = gPars$cols_tr[2], border = NA)
points(nb,
       co[, 11],
       type = 'p',
       pch = 19,
       col = gPars$cols[3])
segments(nb, co[, 12], nb, co[, 13], col = gPars$cols[3])
abline(h = tabResIso$fvX11,
       lty = 2,lwd = 2* gPars$lwd,
       col = gPars$cols[3])
points(nb,
       co[, 14],
       type = 'p',
       pch = 19,
       col = gPars$cols[4])
segments(nb, co[, 15], nb, co[, 16], col = gPars$cols[4])
abline(h = tabResIso$fvX21,
       lty = 2, lwd = 2* gPars$lwd,
       col = gPars$cols[4])
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
intrvOptU  = ErrViewLib::genIntervals(
  uEtrain, nBinOpt,
  equiPop = TRUE,
  logBin  = TRUE)
nBinOptu = intrvOptU$nbr
intrvOptX  = ErrViewLib::genIntervals(
  X1train, nBinOpt,
  equiPop = TRUE,
  logBin  = FALSE)
nBinOptX = intrvOptX$nbr

# Training
savg = sqrt(mean((Etrain/uEtrain)^2))
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

gs = globScore(Etrain, uEtrains, nBinStat, intrvStat, io0, io1, io2)
resp = list()
resp$cal   = gs[1]
resp$cons  = gs[2]
resp$adap1 = gs[3]
resp$adap2 = gs[4]
resp$glob  = sum(gs[1:4])
resp$nll   = gs[5]
resp$ence  = gs[6]
resp$uce   = gs[7]
tabRes2D[['train']] = resp


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

gs = globScore(Etest, uEtests, nBinStat, intrvStat, io0t, io1t, io2t)
resp = list()
resp$cal   = gs[1]
resp$cons  = gs[2]
resp$adap1 = gs[3]
resp$adap2 = gs[4]
resp$glob  = sum(gs[1:4])
resp$nll   = gs[5]
resp$ence  = gs[6]
resp$uce   = gs[7]
tabRes2D[['test']] = resp

sink(file =  file.path(tabDir,paste0('tabRes2D.tex')))
print(knitr::kable(as.data.frame(do.call(rbind,tabRes2D)),'latex'))
sink()
