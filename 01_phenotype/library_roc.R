library(MASS); library(glmpath); library(glmnet)

ROC.Est.FUN <- function(Di, yyi, yy0, fpr0=NULL, wgti=NULL, yes.smooth=F) 
{
  out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
  if(is.null(wgti)){wgti=rep(1, length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));
  mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0
  for(k in 1:pp)
  {
    yy = yy0;
    if(!is.null(fpr0)){
      tpr.all = S.FUN(yyi[,k], Yi=yyi[,k], Di*wgti,yes.smooth=yes.smooth);
      fpr.all = S.FUN(yyi[,k], Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
      TPR = approx(c(0, fpr.all,1), c(0, tpr.all,1),fpr0,method="linear",rule=2)$y;
      TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR);
      yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))
      FPR = S.FUN(yy, Yi=yyi[,k], (1-Di)*wgti, yes.smooth=yes.smooth)
    } else {
      TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth);
      FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
    }
    out.yy = cbind(out.yy, yy)
    out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
    out.TPR = cbind(out.TPR, TPR); out.FPR <- cbind(out.FPR, FPR)
    PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
    out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
    AUC = sum(S.FUN(yyi[,k], Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
    out.AUC <- c(out.AUC, AUC)
  }
  out = c(out.AUC, out.yy, out.pp, out.FPR, out.TPR, out.PPV, out.NPV)
  out
  
}


Sinv.FUN <- function(uu, Yi, Di, yes.smooth=F)
{
  yy0 <- unique(sort(Yi, decreasing = T));
  ss0 <- S.FUN(yy0, Yi, Di, yes.smooth = yes.smooth)
  return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
}


S.FUN <- function(yy,Yi,Di,yes.smooth=F) 
{
  if(yes.smooth){
    Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
    c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  } else {
    return((sum.I(yy,"<", Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  }
}
  

sum.I <- function(yy, FUN, Yi, Vi=NULL)
{
  if (FUN=="<"|FUN==">=") {yy <- -yy; Yi <- -Yi}
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
  if (!is.null(Vi)) {
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else tmpind <- order(Yi)
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}
  

my.roc.fun3 = function(junk, cut.fpr.list=c(0.05,0.1)){
  myauc = junk[1]
  
  tmp = data.frame(matrix(junk[-1], ncol=6, byrow=F))
  colnames(tmp) = c("cut","p.pos", "fpr", "tpr", "ppv","npv")
  tmp=data.frame(sapply(1:ncol(tmp), function(kk) approx(tmp[,"fpr"], tmp[,kk], seq(0.01,0.99,by=0.01), rule=2)$y))
  colnames(tmp) = c("cut","p.pos", "fpr", "tpr", "ppv","npv")
  myroc=NULL
  for(cut.fpr in cut.fpr.list) {
    myroc=c(myroc,as.numeric(tmp[which(abs(round(tmp$fpr,3)-cut.fpr)==min(abs(round(tmp$fpr,3)-cut.fpr)))[length(which(abs(round(tmp$fpr,3)-cut.fpr)==min(abs(round(tmp$fpr,3)-cut.fpr))))],
                                  c("p.pos", "tpr", "ppv","npv")]))
  }
  myroc=c(myauc,myroc)
  names(myroc)=c("auc", unlist(lapply(cut.fpr.list,function(cut.fpr)
    paste0(paste0(c("p.pos", "tpr", "ppv", "npv"),"_fpr"), cut.fpr))))
  myroc
}

VTM<-function(vc,dm) {
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


#ROC.Est.FUN.approx <- function(Di,yyi,wgti=NULL,yes.smooth=F)


