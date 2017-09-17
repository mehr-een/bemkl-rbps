#
# Rule based protein selection
#
# Y: a binary vector of outputs
# X: a binary matrix of features
# 

rm(list=ls())
set.seed(109)
setwd('/Users/meali/Downloads/bemkl-rbps-github/RBPS/')

computeOR <-function(x)
{
  x1 = 0
  for(i in 1:(ncol(x)-1))
    for(j in (i+1):ncol(x)){
      x1 = cbind(x1,x[,i] | x[,j])
      colnames(x1)[ncol(x1)] = paste0(colnames(x)[i]," | ",colnames(x)[j])
    }
  x1 = x1[,-1]
  
  x2 = 0
  for(i in 1:(ncol(x)-1))
    for(j in (i+1):ncol(x)){
      x2 = cbind(x2,x[,i] | !x[,j])
      colnames(x2)[ncol(x2)] = paste0(colnames(x)[i]," | !",colnames(x)[j])
    }
  x2 = x2[,-1]
  
  x3 = 0
  for(i in 1:(ncol(x)-1))
    for(j in (i+1):ncol(x)){
      x3 = cbind(x3,!x[,i] | x[,j])
      colnames(x3)[ncol(x3)] = paste0("!",colnames(x)[i]," | ",colnames(x)[j])
    }
  x3 = x3[,-1]
  
  x4 = 0
  for(i in 1:(ncol(x)-1))
    for(j in (i+1):ncol(x)){
      x4 = cbind(x4,x[,i] | x[,j])
      colnames(x4)[ncol(x4)] = paste0("!",colnames(x)[i]," | !",colnames(x)[j])
    }
  x4 = x4[,-1]
  
  return(list(x=x,xN=!x,x1=x1,x2=x2,x3=x3,x4=x4))
}

ev.sim.real.bidirection <- function(y,x)
{
  inds = which(is.na(y))
  if(length(inds)>0) {y = y[-inds]; x = x[-inds] }
  x[x != 1] = -1
  
  x = as.numeric(x);
  score = 0
  inds = which(x == 1); if(length(inds)!=0) score = mean(y[inds])*(1+ 2*(length(inds)/length(y)) )
  inds = which(x == -1); if(length(inds)!=0) score = score/(mean(y[inds])+1)
  return(score)
}

get.best.feature <- function(dName)
{
  rules = list(x=x, xN = matrix(as.numeric(!x),nrow(x),ncol(x)) )
  nOneSide = 0
  selRules = list()
  selRules$org = array(NA,dim = c(length(dName),ncol(rules$x),2),dimnames = list(dName,colnames(rules$x),c("p","n"))) #p plain, n negated
  for(d in 1:length(dName))
  {
    yd = Y[BC==1,dName[d]]
    if(length(grep(pattern="real",ev.fun.text))>0) {
      maxScore = 0; minScore = Inf;
      for(nu in 1:(length(yd)-1)){
        Ytmp = rep(-1,length(yd)); Ytmp[order(yd,decreasing=T)[1:nu]] = 1 
        tMaxScore <- ev.fun(yd,Ytmp);
        if(tMaxScore > maxScore)
          maxScore <- tMaxScore
        if(tMaxScore < minScore)
          minScore <- tMaxScore
      }
    }
    for(i in 1:ncol(rules$x)){
      selRules$org[d,i,1] = ev.fun(yd,rules$x[,i])*(sum(rules$x[,i] == 1) > nOneSide && sum(rules$x[,i] == 0) > nOneSide)/maxScore
      selRules$org[d,i,2] = ev.fun(yd,rules$xN[,i])*(sum(rules$xN[,i] == 1) > nOneSide && sum(rules$xN[,i] == 0) > nOneSide)/maxScore
    }
  }
  
  rs = selRules$org[1,,]
  rownames(rs) = gsub(rownames(rs),pattern="RPPA",replacement = "")
  rownames(rs) = gsub(rownames(rs),pattern="MS",replacement = "")
  rind = which.max(apply(rs,1,max))
  cind = which.max(apply(rs,2,max))
  if(names(cind)=="p")
    nm = paste0(rownames(rs)[rind],"_up")
  if(names(cind)=="n")
    nm = paste0(rownames(rs)[rind],"_down")
  
  b.feat = rs[rind,cind]
  names(b.feat) = nm
  return(b.feat)
}

get.top.rules <- function(dNames, nTop)
{
  nOneSide = 2
  selRules = list()
  selRules$org = array(NA,dim = c(length(dNames),ncol(rules$x),2),dimnames = list(dNames,colnames(rules$x),c("p","n"))) 
  selRules$int = array(NA,dim = c(length(dNames),ncol(rules$x1),4),dimnames = list(dNames,colnames(rules$x1),c("p","1","2","n")))
  for(d in 1:length(dNames))
  {
    yd = Y[,dNames[d]]
    if(length(grep(pattern="real",ev.fun.text))>0) {
      maxScore = 0; minScore = Inf;
      for(nu in 2:(length(yd)-2)){
        Ytmp = rep(-1,length(yd)); Ytmp[order(yd,decreasing=T)[1:nu]] = 1 
        tMaxScore <- ev.fun(yd,Ytmp);
        if(tMaxScore > maxScore)
          maxScore <- tMaxScore
        if(tMaxScore < minScore)
          minScore <- tMaxScore
      }
    }
    
    for(i in 1:ncol(rules$x1)){
      selRules$int[d,i,1] = ev.fun(yd,rules$x1[,i])*(sum(rules$x1[,i] == 1) > nOneSide && sum(rules$x1[,i] == 0) > nOneSide)/maxScore
      selRules$int[d,i,2] = ev.fun(yd,rules$x2[,i])*(sum(rules$x2[,i] == 1) > nOneSide && sum(rules$x2[,i] == 0) > nOneSide)/maxScore
      selRules$int[d,i,3] = ev.fun(yd,rules$x3[,i])*(sum(rules$x3[,i] == 1) > nOneSide && sum(rules$x3[,i] == 0) > nOneSide)/maxScore
      selRules$int[d,i,4] = ev.fun(yd,rules$x4[,i])*(sum(rules$x4[,i] == 1) > nOneSide && sum(rules$x4[,i] == 0) > nOneSide)/maxScore
    }
    for(i in 1:ncol(rules$x)){
      selRules$org[d,i,1] = ev.fun(yd,rules$x[,i])*(sum(rules$x[,i] == 1) > nOneSide && sum(rules$x[,i] == 0) > nOneSide)/maxScore
      selRules$org[d,i,2] = ev.fun(yd,rules$xN[,i])*(sum(rules$xN[,i] == 1) > nOneSide && sum(rules$xN[,i] == 0) > nOneSide)/maxScore
    }
  }
  colnames(selRules$org)[grep(colnames(selRules$org),pattern = "FH")]
  
  tt = nTop
  f.rules = list()
  for(d in 1:length(dNames))
  {
    f.rules[[d]] = rep(NA,tt)
    tinds = order(apply(round(selRules$int[d,,],2),1,max),decreasing=T)[1:tt]
    for(it in 1:length(tinds)){
      tp = names(which.max(selRules$int[d,tinds[it],]))
      nm = rownames(selRules$int[d,tinds,])[it]
      nm = unlist(strsplit(nm,split=" | ",fixed=T))
      nm = gsub(nm,pattern="MS",replacement = "")
      nm = gsub(nm,pattern="RPPA",replacement = "")
      sc = round(max(selRules$int[d,tinds[it],]),2)
      
      if(tp == "p"){
        nm[1] = paste0(nm[1],"_up")
        nm[2] = paste0(nm[2],"_up")
      }
      if(tp == "1"){
        nm[1] = paste0(nm[1],"_up")
        nm[2] = paste0(nm[2],"_down")
      }
      if(tp == "2"){
        nm[1] = paste0(nm[1],"_down")
        nm[2] = paste0(nm[2],"_up")
      }
      if(tp == "n"){
        nm[1] = paste0(nm[1],"_down")
        nm[2] = paste0(nm[2],"_down")
      }
      nm = paste0(nm,collapse=" OR ")
      f.rules[[d]][it] = sc
      names(f.rules[[d]])[it] = nm
    }
  }
  names(f.rules) = dNames
  
  return(f.rules)
}


#### load data
load("ProtViews.RData")
load("DrugResponseTargeted.RData")

rownames(targeted_response) = targeted_agents
Y = t(targeted_response)

#select X 
rownames(prot_bin) = colnames(prot)
x = prot_bin
x[x == -1] = 0

#### run analysis
rules = computeOR(x)
ev.fun.text = "ev.sim.real.bidirection"
ev.fun = get(ev.fun.text)

dNames = c("Alvespimycin","Macbecin II","Selumetinib","Tamoxifen","Tanespimycin","Alvocidib")
top.rules = get.top.rules(dNames, nTop=10)

print("Features for Alvespimycin")
print(top.rules$Alvespimycin)
# HSP90AA1 OR HSP90AB1 are primary targets

print("Features for Macbecin II")
print(top.rules$'Macbecin II')
# HSP90AA1 OR HSP90AB1 are primary targets

print("Features for Selumetinib")
print(top.rules$Selumetinib)
# No primary target. Though NPM1 feature combinations occur most often, hence selecting the NPM1 feature with highest score.

print("Features for Tamoxifen")
print(top.rules$Tamoxifen)
# No primary target. Though DDX5 feature combinations occur most often, hence selecting the DDX5 feature with highest score.

print("Features for Tanespimycin")
print(top.rules$Tanespimycin)
# HSP90AA1 OR HSP90AB1 are primary targets

print("Features for Alvocidib")
print(top.rules$Alvocidib)
# No primary target. Though NPM1 feature combinations occur most often, hence selecting the NPM1 feature with highest score.

#### run lapatinib
ctype = read.table("CellTypes.csv",sep=";",header=F)
ct = ctype[,2]

BC = as.numeric(ct == "Breast" | ct == "Ovarian");
x = x[BC==1,] 
x = x[,apply(x,2,sd)>0]

dName = "Lapatinib"
b.feat = get.best.feature(dName)
print("Features for Lapatinib")
print(b.feat)
# Best feature for Lapatinib

