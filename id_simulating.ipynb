import numpy as np
import dynamo as dyn
import matplotlib.pyplot as plt
import skdim
import time


simulator = dyn.sim.BifurcationTwoGenes(dyn.sim.bifur2genes_params, tau=5)
simulator.simulate([0,50], n_cells=5000)
adata = simulator.generate_anndata()

adata.obsm['X_raw'] = adata.layers['X']
adata.obsm['velocity_raw'] = adata.layers['V']

# dyn.pl.zscatter(adata, basis='raw', color='trajectory', cmap='Set1')
# dyn.pl.zstreamline(adata, basis='raw')
# plt.axis('on')
# plt.show()

start = time.time()
lpca = skdim.id.lPCA()
lpca_id = lpca.fit_transform(adata.obs)
end = time.time()
print("lPCA: " + str(lpca_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
mind_ml = skdim.id.MiND_ML()
mind_ml_id = mind_ml.fit_transform(adata.obs)
end = time.time()
print("MiND_ML: " + str(mind_ml_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
mada = skdim.id.MADA()
mada_id = mada.fit_transform(adata.obs)
end = time.time()
print("MADA: " + str(mada_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
ess = skdim.id.ESS()
ess_id = ess.fit_transform(adata.obs)
end = time.time()
print("ESS: " + str(ess_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
knn = skdim.id.KNN()
knn_id = knn.fit_transform(adata.obs)
end = time.time()
print("KNN: " + str(knn_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
corr = skdim.id.CorrInt()
corr_id = corr.fit_transform(adata.obs)
end = time.time()
print("CorrInt: " + str(corr_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
danco = skdim.id.DANCo()
danco_id = danco.fit_transform(adata.obs)
end = time.time()
print("DanCo: " + str(danco_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
fisher = skdim.id.FisherS()
fisher_id = fisher.fit_transform(adata.obs)
end = time.time()
print("FischerS: " + str(fisher_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
mle = skdim.id.MLE()
mle_id = mle.fit_transform(adata.obs)
end = time.time()
print("MLE: " + str(mle_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
mom = skdim.id.MOM()
mom_id = mom.fit_transform(adata.obs)
end = time.time()
print("MOM: " + str(mom_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
tle = skdim.id.TLE()
tle_id = tle.fit_transform(adata.obs)
end = time.time()
print("TLE: " + str(tle_id) + "\tin {:.2f}s".format(end-start))

start = time.time()
twoNN = skdim.id.TwoNN()
twoNN_id = twoNN.fit_transform(adata.obs)
end = time.time()
print("twoNN: " + str(twoNN_id) + "\tin {:.2f}s".format(end-start))

