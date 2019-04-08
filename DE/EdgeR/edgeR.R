library(edgeR)
library('Matrix')
setwd('/data/yosef2/users/chenling/CSF')
# celltypes 'B1', 'B2', 'CD4', 'CD8a', 'CD8n', 'Gran', 'MegaK', 'Mono', 'NK1', 'NK2', 'Tdg', 'Tregs', 'mDC1', 'mDC2', 'ncMono', 'pDC', 'plasma'
i = commandArgs(trailingOnly=T)[1]

count = readMM('RFiles/all_data.mtx')
count = t(count)

genenames = read.csv('RFiles/genenames.csv',header=F)[,1]
genenames = as.character(genenames)


batchid = read.csv('RFiles/batchid.csv',header=F)[,1]

isMS = read.csv('RFiles/isMS.csv',header=F)[,1]
isCSF = read.csv('RFiles/isCSF.csv',header=F)[,1]
state = rep('MS',length(isMS))
state[isMS=='False']='control'
tissue = rep('CSF',length(isCSF))
tissue[isCSF=='False']='PBMC'
state = factor(state)
tissue = factor(tissue)
celllabels = read.csv('RFiles/celllabels.csv',header=F)[,1]
# celllabels = factor(celllabels)


y <- DGEList(counts=count, group=(celllabels==i))
y$genes <- data.frame(Symbol=genenames)
cdr <- scale(colMeans(count > 0))
istype = as.factor(celllabels==i) 
design= model.matrix(~ cdr  + as.factor(batchid) + istype )
rank = qr(design)$pivot[seq_len(qr(design)$rank)]
design=  design[, rank]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
# contrast = rep(-1/(nclasses-1), nclasses)
# contrast[i] = 1
# qlf <- glmQLFTest(fit, contrast=contrast)
qlf <- glmQLFTest(fit)
write.csv(qlf$table, file=paste('EdgeR/allcluster.batchcorrected.',i,'.edgeR.csv',sep=''))
