import scipy.io
from scipy.cluster import hierarchy as hc
from scipy.spatial import distance as scd
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np
import fastcluster as fc
import csv

def diss_cells(M,j):
	print 'Computing dissimilarity matrix for cells from correlation matrix'
	c = pairwise_distances(X=M.T.toarray(),metric='correlation',n_jobs=j)
	np.fill_diagonal(c,0.)
	c = scd.squareform(c,force='tovector',checks=False)
	c = np.nan_to_num(c)
	return c

def diss_genes(M,j):
	print 'Computing dissimilarity matrix for genes from correlation matrix'
	c = pairwise_distances(X=M.toarray(),metric='correlation',n_jobs=j)
	np.fill_diagonal(c,0.)
	c[np.isnan(c)]=1
	c=scd.squareform(c,force='tovector',checks=False)
	return c

def create_linkage(diss,method):
	print 'Computing linkage'
	z = fc.linkage(diss,method=method)
	return z.clip(min=0)

def create_leaves_list(z):
	return hc.leaves_list(z)

def save_matrix(M,output):
	scipy.io.mmwrite(output,M)


## read files
M = scipy.io.mmread('matrix.mtx').tocsc()
labels_cells = [row[0] for row in csv.reader(open('barcodes.tsv'), delimiter="\t")]
labels_genes = [row[1].upper() for row in csv.reader(open('genes.tsv'), delimiter="\t")]

## clustering cells (columns)
diss_cells = diss_cells(M,-1)
z_cells = create_linkage(diss_cells,'single')
idx_cells = create_leaves_list(z_cells)

## clustering genes (rows)
diss_genes = diss_genes(M,-1)
z_genes = create_linkage(diss_genes,'single')
idx_genes = create_leaves_list(z_genes)

## sorting matrix and labels
M = M[idx_genes,:]
M = M[:,idx_cells]
hclust_cells = [labels_cells[int(i)] for i in idx_cells]
hclust_genes = [labels_genes[int(i)] for i in idx_genes]

## plotting heatmap
#bash command: pip install pydendroheatmap
import pydendroheatmap as pdh
heatmap = pdh.DendroHeatMap(heat_map_data=M.toarray(), left_dendrogram=z_genes, top_dendrogram=z_cells)
heatmap.colormap = heatmap.yellowBlackBlue
heatmap.row_labels = hclust_genes
heatmap.show()

## Matlab can plot large heatmaps better than any other program.
## I strongly advise to load the clustered matrix in Matlab
## and use imagesc() to visualize the results. Students can download
## Matlab for free from the access.caltech.edu website.
save_matrix(M,'clustered.mtx')










