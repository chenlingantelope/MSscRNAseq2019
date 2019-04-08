
library('Matrix')
library(devtools)
load_all("VISION")
library(RcppCNPy)

doubletscore <- npyLoad("RFiles/doublet_score.npy")
merged_exprs = readMM('RFiles/CD4data.mtx')
merged_exprs = t(merged_exprs)
exprs = as.matrix(merged_exprs)

batchid = read.csv('RFiles/batchid.csv',header=F)[,1]
isMS = read.csv('RFiles/isMS.csv',header=F)[,1]
isCSF = read.csv('RFiles/isCSF.csv',header=F)[,1]
clusters = read.csv('RFiles/celllabels.csv',header=F)[,1]
CD4clusters = read.csv('RFiles/CD4labels.csv',header=F)[,1]
cellsize = colSums(exprs)
ngenes = colSums(exprs>0)

meta = data.frame(
	batch = as.factor(batchid)[clusters=='CD4'], 
	isMS = as.factor(isMS)[clusters=='CD4'], 
	isCSF = as.factor(isCSF)[clusters=='CD4'],
	clusters= as.factor(CD4clusters), 
	cellsize = cellsize, 
	ngenes = ngenes, 
	doubletscore = doubletscore[clusters=='CD4'])

colnames(meta) = c('batch','isMS','isCSF','clusters','cellsize','ngenes','doubletscore')
genenames = read.csv('RFiles/genenames.csv',header=F,as.is=T)[,1]

latent_u = read.csv('RFiles/latent_u.csv',header=F)
latent = read.csv('RFiles/latent.csv',header=F)

latent_u = latent_u[clusters=='CD4',]
latent = latent[clusters=='CD4',]

raw_exprs = exprs
n.umi = median(colSums(exprs))
exprs = apply(exprs, 2, function(x) (x * n.umi) / sum(x))

rownames(exprs) = genenames
colnames(exprs) = as.character(c(1:length(exprs[1,])))
rownames(raw_exprs) = genenames
colnames(raw_exprs) = as.character(c(1:length(exprs[1,])))
rownames(latent) = as.character(c(1:length(exprs[1,])))
rownames(latent_u) = as.character(c(1:length(exprs[1,])))

gmt = c('Hemato.geneset.gmt','cell_cycle_Tirosh.gmt','signatures_NY_private.gmt','c7.all.v6.2.symbols.gmt')
gmt = sapply(gmt,function(x){paste('signatures/',x,sep='')})
# -- find genes to use for projections
# f.genes = VISION:::filterGenesFano(exprs)
# -- Create a VISION object with existing matrix

vis <- Vision(data = exprs + 1, nomodel = TRUE,
                  unnormalizedData = raw_exprs,
                  signatures = gmt,
                  meta = meta,
                  pool = F,
                  latentSpace = latent,
                  cellsPerPartition = 10)

vis <- addProjection(vis, "UMAP", latent_u)
vis <- analyze(vis)

saveRDS(vis,file = "Vision.CD4.rds")

# vis = readRDS("Vision.withmeta.allcells.rds")
# viewResults(vis,host='0.0.0.0',browser=F,port=8111)

# pooledvis = readRDS('Vision.simple.10pool.rds')
# viewResults(pooledvis,host='0.0.0.0',browser=F,port=8112)


