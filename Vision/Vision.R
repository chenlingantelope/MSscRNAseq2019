
library('Matrix')
library(devtools)
load_all("VISION")
library(RcppCNPy)

doubletscore <- npyLoad("Notebooks/meta/doublet_score.npy")
merged_exprs = readMM('RFiles/all_data.mtx')
merged_exprs = t(merged_exprs)
exprs = as.matrix(merged_exprs)

batchid = read.csv('RFiles/batchid.csv',header=F)[,1]
isMS = read.csv('RFiles/isMS.csv',header=F)[,1]
isCSF = read.csv('RFiles/isCSF.csv',header=F)[,1]
clusters = read.csv('RFiles/celllabels.csv',header=F,sep=',')[,1]
cellsize = colSums(exprs)
ngenes = colSums(exprs>0)

meta = data.frame(batch = as.factor(batchid), isMS = as.factor(isMS), isCSF = as.factor(isCSF),
	clusters= as.factor(clusters), cellsize = cellsize, ngenes = ngenes, doubletscore = doubletscore)

colnames(meta) = c('batch','isMS','isCSF','clusters','cellsize','ngenes','doubletscore')
genenames = read.csv('RFiles/genenames.csv',header=F,as.is=T)[,1]

latent_u = read.table('RFiles/latent_u.csv',sep=',')
latent = read.table('RFiles/latent.csv',sep=',')

raw_exprs = exprs
n.umi = median(colSums(exprs))
exprs = apply(exprs, 2, function(x) (x * n.umi) / sum(x))

rownames(exprs) = genenames
colnames(exprs) = as.character(c(1:length(isMS)))
rownames(raw_exprs) = genenames
colnames(raw_exprs) = as.character(c(1:length(isMS)))
rownames(latent) = as.character(c(1:length(isMS)))
rownames(latent_u) = as.character(c(1:length(isMS)))

gmt = list.files('signatures/')
gmt = sapply(gmt, function(X){paste('signatures/',X,sep='')})
# -- find genes to use for projections
# f.genes = VISION:::filterGenesFano(exprs)
# -- Create a VISION object with existing matrix
vis <- Vision(data = exprs + 1, nomodel = TRUE,
                  unnormalizedData = raw_exprs,
                  signatures = gmt,
                  meta = meta,
                  pool = F,
                  latentSpace = latent,
                  cellsPerPartition = 1)

vis <- addProjection(vis, "UMAP", latent_u)

options(mc.cores=10)
vis <- analyze(vis)

saveRDS(vis,file = "VisionRDS/Vision.rds")

viewResults(vis,host='0.0.0.0',browser=F,port=8111)

# pooledvis = readRDS('Vision.simple.10pool.rds')
# viewResults(pooledvis,host='0.0.0.0',browser=F,port=8112)
