import scanpy as sc
import scvi
from scvi.model.utils import mde

adata = sc.read_h5ad('data/allcell_seu.h5ad')

adata.raw = adata  # keep full dimension safe
sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="species",
    subset=True,
)

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="species")
vae = scvi.model.SCVI(adata, gene_likelihood="nb", n_layers=3, n_latent=30)

vae.train()
adata.obsm["scVI"] = vae.get_latent_representation()

sc.pp.neighbors(adata, use_rep="scVI")
adata.obsm["X_mde"] = mde(adata.obsm["scVI"])

vae.save('model/')
adata.write('tmp/allcell.h5ad')
