for f in {0..11}
do
	Rscript EdgeR/edgeR.CD4.R  $f
done


for f in 'plasma','B1','B2','CD4','Tregs','Tdg','CD8n','CD8a','NK1','NK2','Mono','ncMono','Gran','mDC1','mDC2','pDC':
do
	do Rscript EdgeR/edgeR.R  $f
done

