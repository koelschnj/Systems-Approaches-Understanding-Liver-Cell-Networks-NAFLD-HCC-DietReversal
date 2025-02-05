import liana as li
import scanpy as sc
import anndata as ad
import pandas as pd
import plotnine as p9
from plotnine import ggplot, aes, geom_point, facet_grid, labs, guide_legend, guides, theme, element_text, element_line, element_rect, theme_set, theme_void

adata = ad.read_h5ad("/lustre/home/koelschnj/RDn2024.h5ad")

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.raw.X
adata

li.mt.show_methods()

# import liana's rank_aggregate
from liana.mt import rank_aggregate

rank_aggregate.describe()

from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean

li.rs.show_resources

li.mt.rank_aggregate(adata,
                     groupby='Abbreviated',
                     resource_name='mouseconsensus',
                     expr_prop=0.1, use_raw=False,
                     verbose=True)

adata.uns['liana_res'].head()

rank_aggregate.describe()

adata.uns['liana_res'].to_csv('LIANA_res_RDn2024.csv', index=False)

li.pl.dotplot(adata = adata,
              colour='magnitude_rank',
              size='specificity_rank',
              inverse_size=True,
              inverse_colour=True,
              top_n=10,
              orderby='magnitude_rank',
              orderby_ascending=True,
              figure_size=(8, 7)
             )

my_plot = li.pl.dotplot(adata = adata,
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
			top_n=10,
			orderby='magnitude_rank',
             		orderby_ascending=True,
			source_labels=["Mono"],
			target_labels=["Mono","Myofibro","Fibro","Hep","Cancer","Stromal"],
              		figure_size=(6, 4),
                        filter_fun=lambda x: x['specificity_rank'] <= 0.01,
                       )
my_plot

my_plot.save("dotplotfilteredRDn2024.png")

adj_my_plot = (my_plot +
 # change theme
 p9.theme_dark() +
 # modify theme
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11,color="#222222"),
     axis_text_x=element_text(color="#222222"),
     axis_text_y=element_text(color="#222222"),
     figure_size=(6, 4)) +
 labs(title="Source")
)

adj_my_plot.save("p9dotplotfilteredRDn2024.png")

my_plot1 = li.pl.dotplot(adata = adata,
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
			orderby='magnitude_rank',
             		orderby_ascending=True,
			source_labels=["Mono"],
			target_labels=["Mono","Myofibro","Fibro","Hep","Cancer","Stromal"],
              		figure_size=(6, 4),
			ligand_complex='Fn1',
			receptor_complex = ['Cd44','Itgav_Itgb1','Sdc4','Itga4_Itgb1'],
                       )
my_plot1

my_plot1.save("dotplotFn1RDn2024.png")

adj_my_plot1 = (my_plot1 +
 # change theme
 p9.theme_dark() +
 # modify theme
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11,color="#222222"),
     axis_text_x=element_text(color="#222222"),
     axis_text_y=element_text(color="#222222"),
     figure_size=(6, 4)) +
 labs(title="Source")
)

adj_my_plot1.save("p9Fn1dotplotfilteredRDn2024.png")

my_plot2 = li.pl.dotplot(adata = adata,
                        colour='magnitude_rank',
                        inverse_colour=True,
                        size='specificity_rank',
                        inverse_size=True,
			orderby='magnitude_rank',
             		orderby_ascending=True,
			source_labels=["Mono", "Cancer","Stromal"],
			target_labels=["Mono","Myofibro","Fibro","Hep","Cancer","Stromal"],
              		figure_size=(13, 5),
			ligand_complex=['F2','Plg','Pros1','Vtn','Nrg4','Angptl4','Hgf'],
			receptor_complex=['Pard3','Mertk','Itgav_Itgb1','Itgav_Itgb5','Met','Sdc2','Sdc4'],
                       )
my_plot2


adj_my_plot2 = (my_plot2 +
 # change theme
 p9.theme_dark() +
 # modify theme
 p9.theme(
     # adjust facet size
     strip_text=p9.element_text(size=11,color="#222222"),
     axis_text_x=element_text(color="#222222"),
     axis_text_y=element_text(color="#222222"),
     figure_size=(13, 5)) +
 labs(title="Source")
)

adj_my_plot2.save("p9S10dotplotfilteredRDn2024.png")



