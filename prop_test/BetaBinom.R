library('VGAM')
# library('rjags')
# library('overlapping')
library('aod')
options(scipen=999)

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
count = read.table(paste(filename, '.count.txt',sep=''))
N = as.numeric(colSums(count))
tissue=c(rep('CSF',8),rep('PBMC',10))
states = c('MS','MS','MS','MS','control','control','control','control','MS','MS','MS','MS','MS','control','control','control','control','control')

# BetaBinomial1 = "
# model{
#     for(i in 1:N){
#     	x[i] ~ dbinom(p1[i]*c[i]+p0[i]*(1-c[i]),n[i])	
# 	    p0[i] ~ dbeta(a0,b0) T(,1-1e-100)
#     	p1[i] ~ dbeta(a1,b1) T(,1-1e-100)
#     }
#     a0 ~ dgamma(1e-10, 1e-10)
#     b0 ~ dgamma(1e-10, 1e-10)
#     a1 ~ dgamma(1e-10, 1e-10)
#     b1 ~ dgamma(1e-10, 1e-10)
# }"


BetaBinomialReg = function(count,N,C,filename){
	res = lapply(c(1:length(count[,1])),function(i){
		x = as.numeric(count[i,])
		# fit = vglm(t(rbind( x, N-x )) ~ C, betabinomialff)
		# p = coef(summary(fit))[c(3,4),4]
		# pred = fitted(fit)
		# fc = pred[C==T][1] / pred[C==F][1]
		data = data.frame(x=x,N=N,C=as.factor(C))
		fm1 <- betabin(cbind(x, N - x) ~ C, ~ 1, data)
		res = summary(fm1)@Coef
		p = res[2,4]
		# pred = fitted(fm1)
		# fc = pred[C==T][1] / pred[C==F][1]
		fc = mean(x[C==T]/N[C==T])/mean(x[C==F]/N[C==F]) 
		return(c(p,fc))
	})
	res = do.call(cbind,res)
	write.table(t(res),col.names=F,row.names=F,quote=F,sep=',',file=filename)	
	return(res)
}



# BinomialReg = function(count,N,C,filename){
# 	res = lapply(c(1:length(count[,1])),function(i){
# 		x = as.numeric(count[i,])
# 		df <- data.frame(
# 		    TotalCells = N,
# 		    numInCluster = x,
# 		    Genotype = as.factor(C)
# 		)
# 		result <- glm(numInCluster/TotalCells~Genotype,
# 		              data=df, family=binomial,
# 		              weights=TotalCells)
# 		p = summary(result)$coef[2,4]
# 		return(c(p,mean(x[C==T]/N[C==T])/mean(x[C==F]/N[C==F]) ))
# 	})
# 	res = do.call(cbind,res)
# 	write.table(t(res),col.names=F,row.names=F,quote=F,sep=',',file=filename)	
# 	return(res)
# }

# ComputeDP = function(count,N,C,filename){
# 	pdf(paste(filename,'.pdf',sep=''))
# 	OV = lapply(c(1:length(count[,1])),function(i){
# 		X = count[i,]
# 		model <- jags.model(textConnection(BetaBinomial1),
# 		                    data = list('x' = as.numeric(X),
# 		                    	'n' = as.numeric(N),
# 		                    	'c' = C,
# 		                    	'N' = length(N)))
# 		update(model, 5000)
# 		params <- jags.samples(model,c('a0','a1','b0','b1'),50000)
# 		mu0 = params[[1]][1,,1]/(params[[1]][1,,1]+params[[3]][1,,1])
# 		mu1 = params[[2]][1,,1]/(params[[2]][1,,1]+params[[4]][1,,1])
# 		par(mfrow=c(3,3))
# 		lower = min(c(mu0,mu1))
# 		upper = max(c(mu0,mu1))
# 		step = (upper-lower)/100
# 		hist(mu0,breaks=seq(lower,upper,step),col='blue',border='blue')
# 		hist(mu1,breaks=seq(lower,upper,step),col='red',add=T,border='red')
# 		res = c(overlap(list(log(mu0),log(mu1)))$OV, mean(mu1)/mean(mu0))
# 		return(list(res,mu0,mu1))
# 	})
# 	mu0 = lapply(OV,function(X){X[[2]]})
# 	mu1 = lapply(OV,function(X){X[[3]]})
# 	OV = sapply(OV,function(X){X[[1]]})
# 	plot(y=-log(OV[1,]+1e-10,10),x=log(OV[2,],10))
# 	dev.off()
# 	write.table(t(OV),col.names=F,row.names=F,quote=F,sep=',',file=filename)
# 	return(list(OV,mu0,mu1))
# }

# OV1 = ComputeDP(count,N,tissue=='CSF','Bayesian.CSF.csv')
# OV2 = ComputeDP(count,N,states=='MS','Bayesian.MS.csv')
# OV3 = ComputeDP(count[,tissue=='CSF'],N[tissue=='CSF'],states[tissue=='CSF']=='MS','Bayesian.CSFMS.csv')
# OV4 = ComputeDP(count[,tissue=='PBMC'],N[tissue=='PBMC'],states[tissue=='PBMC']=='MS','Bayesian.PBMCMS.csv')
# OV5 = ComputeDP(count[,states=='control'],N[states=='control'],tissue[states=='control']=='CSF','Bayesian.CSF_control.csv')

# tres = c(4.29382293, 3.33025022, 7.61345436, 0.97755275, 0.91960306,
#        1.5547467 , 0.36301796, 0.78215611, 0.57243582, 0.46252548,
#        5.18488324, 2.8403952 , 0.78215611, 0.59855511, 2.02233809,
#        0.91651382, 0.46077512, 1.2522841 , 2.47588862, 0.21768479,
#        0.71656765, 8.59107891)

# count[12, states=='MS' & tissue=='CSF']
# res1 = BinomialReg(count,N,tissue=='CSF','count_test/binomreg.CSF.csv')
# res2 = BinomialReg(count,N,states=='MS','count_test/binomreg.MS.csv')
# res3 = BinomialReg(count[,tissue=='CSF'],N[tissue=='CSF'],states[tissue=='CSF']=='MS','count_test/binomreg.CSFMS.csv')
# res4 = BinomialReg(count[,tissue=='PBMC'],N[tissue=='PBMC'],states[tissue=='PBMC']=='MS','count_test/binomreg.PBMCMS.csv')
# res5 = BinomialReg(count[,states=='control'],N[states=='control'],tissue[states=='control']=='CSF','count_test/binomreg.CSF_control.csv')


res1 = BetaBinomialReg(count,N,tissue=='CSF',paste(filename,'.betabinomreg.CSF.csv',sep=''))
res2 = BetaBinomialReg(count,N,states=='MS',paste(filename,'.betabinomreg.MS.csv',sep=''))
res3 = BetaBinomialReg(count[,tissue=='CSF'],N[tissue=='CSF'],states[tissue=='CSF']=='MS',
	paste(filename,'.betabinomreg.CSFMS.csv',sep=''))
res4 = BetaBinomialReg(count[,tissue=='PBMC'],N[tissue=='PBMC'],states[tissue=='PBMC']=='MS',
	paste(filename,'.betabinomreg.PBMCMS.csv',sep=''))
res5 = BetaBinomialReg(count[,states=='control'],N[states=='control'],tissue[states=='control']=='CSF',
	paste(filename,'.betabinomreg.CSF_control.csv',sep=''))
