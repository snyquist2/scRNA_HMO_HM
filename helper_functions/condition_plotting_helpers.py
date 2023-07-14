import scanpy as sc
import anndata
import matplotlib.pylab as plt
import sys
import pandas as pd
#sys.path.insert(0,'/data/cb/nyquist/breast_milk/breastMilk/')
import plotting_helpers as hh
import random
import numpy as np
import scipy.sparse
import seaborn as sns

def split_umap_by_category(adata, groupby, colorby=None,gene=None, nrows=None, ncols=None, figsize_multiplier=4, markersize=.3, color_dict = None,other_color="lightgray",gene_cmap="Reds", dim_reduction="umap",cat_order=None):
    '''
    plot one umap by each category in groupby and color the plot by categories in color_by

    if nrows or ncols is none, set them to 1 row and num unique categories cols

    color_dict will be set to the adata color parameter if it is none
    '''


    if nrows is None or ncols is None:
        nrows = 1
        ncols = len(adata.obs[groupby].unique())
    if colorby is None:
        colorby = groupby
    

    if color_dict is None and gene is None:
        color_dict = dict(zip(adata.obs[colorby].cat.categories,adata.uns[colorby+"_colors"]))
    if gene is not None:
        if scipy.sparse.issparse(adata.raw[:,gene].X):
            genedf = pd.DataFrame(adata.raw[:,gene].X.todense(),index=adata.obs_names)
        else:
            genedf = pd.DataFrame(adata.raw[:,gene].X,index=adata.obs_names)
        gene_min = genedf.min()
        gene_max = genedf.max()
    if gene is None: 
        fig,ax = plt.subplots(nrows,ncols,   figsize=(figsize_multiplier*ncols,figsize_multiplier*nrows))
    else:
        fig,ax = plt.subplots(nrows,ncols,   figsize=(figsize_multiplier*(ncols+.5),figsize_multiplier*nrows))
    if cat_order is not None:
        vals=cat_order
    else:
        vals = sorted(list(adata.obs[groupby].unique()))
    val_to_coords = {}
    val_to_exp={}
    for ind,val in enumerate(vals):
        val_to_coords[val] = adata.obsm["X_"+dim_reduction][adata.obs[groupby]==val]
        if gene is not None:
            val_to_exp[val] = genedf.loc[adata.obs[groupby]==val].values

    for ind,val in enumerate(vals):
        other_vals = list(vals)
        other_vals.remove(val)
        row = int(ind/ncols)
        col = ind%ncols
        if nrows >1:
            this_ax = ax[row,col]
        else:
            this_ax = ax[col]
        
        for v in other_vals:
            this_ax.scatter(val_to_coords[v][:,0],val_to_coords[v][:,1],c=other_color,s=markersize,marker="o",alpha=1)
        
        adata_sub = adata[adata.obs[groupby]==val]
        vals_2 = list(adata_sub.obs[colorby].unique())
        val_to_coords2={}
        for ind2,val_2 in enumerate(vals_2):
            val_to_coords2[val_2] = adata_sub.obsm["X_"+dim_reduction][adata_sub.obs[colorby]==val_2]
            if gene is None:
                this_ax.scatter(val_to_coords2[val_2][:,0],val_to_coords2[val_2][:,1],s=markersize,marker="o",c =color_dict[val_2],alpha=1 )
            else:
                m=this_ax.scatter(val_to_coords2[val_2][:,0],val_to_coords2[val_2][:,1],s=markersize,marker=".",c = val_to_exp[val_2].T.flatten(),cmap=gene_cmap,vmin=gene_min,vmax=gene_max)
                if ind == len(adata.obs[groupby].unique())-1:
                    plt.colorbar(m,ax=ax)
        this_ax.set_title(val)
        this_ax.set_xticks([])
        this_ax.set_yticks([])
    if gene is not None:
        fig.suptitle(gene)
    return fig


def make_downsampled_heatmap(adata, cells_per_group, cluster_key, n_genes=5, key=None):
    sc.tl.dendrogram(adata,groupby=cluster_key)
    sampled_cellnames = []
    for i in adata.obs[cluster_key].unique():
        poss_cells = list(adata[adata.obs[cluster_key]==i].obs_names)
        sampled_cellnames += random.sample(poss_cells, min(cells_per_group, len(poss_cells)))
    if key is None:
        key="rank_genes_groups"
    sc.pl.rank_genes_groups_heatmap(adata[sampled_cellnames], groupby=cluster_key, n_genes=n_genes, swap_axes=True, show_gene_labels=True,cmap=cmap, standard_scale="var",show=False,key=key)



def grouped_stacked_bars(adata, plot_sep, x_value, color_value,normalized=True, figwidth_mult = 6, figheight=7, plot_sep_cat_order=[],color_order=None):
    nplots = len(adata.obs[plot_sep].unique())
    plot_no = 0
    adata.obs[plot_sep] = adata.obs[plot_sep].astype("category")
    if len(plot_sep_cat_order) ==  0:
        plot_sep_cat_order = adata.obs[plot_sep].cat.categories
    fig,ax = plt.subplots(1,nplots,  sharey='all', figsize=(figwidth_mult*nplots,figheight),gridspec_kw={'width_ratios': [len(set(adata.obs.loc[adata.obs[plot_sep] == p,x_value])) for p in plot_sep_cat_order]})
    for p in plot_sep_cat_order:
        a = adata[adata.obs[plot_sep]==p]
        if normalized:
            hh.normalized_stacked_bar_plot(a, x_value, color_value,  legend=False, ax=ax[plot_no], color_order=color_order)
        else:
            hh.stacked_bar_plot(a, x_value, color_value,  legend=False, ax=ax[plot_no])
        ax[plot_no].set_title(p)
        plot_no += 1
    name_to_color={}
    for a in ax:
        handles, labels = a.get_legend_handles_labels()
        if type(color_order) != type(None):
            for h in color_order:
                name_to_color[h] = handles[labels.index(h)]
        for i,h in enumerate(labels):
            name_to_color[h] = handles[i]
        a.tick_params(left=False)
        a.spines['right'].set_visible(False)
        a.spines['top'].set_visible(False)
        a.spines['left'].set_visible(False)
    fig.legend(reversed(list(name_to_color.values())), reversed(list(name_to_color.keys())), loc='center right', ncol=1,bbox_to_anchor=(1.1, .7),borderaxespad=0)
    plt.subplots_adjust(wspace=.1, bottom=.3)

def volcano_plot(adata,col,pmin,lmin,lmax):
    df=pd.concat([pd.DataFrame(adata.uns["rank_genes_groups"]["names"])[col],pd.DataFrame(adata.uns["rank_genes_groups"]["logfoldchanges"])[col],pd.DataFrame(adata.uns["rank_genes_groups"]["pvals_adj"])[col],pd.DataFrame(adata.uns["rank_genes_groups_filtered"]["names"])[col]],axis=1)
    df.columns=["names","logfoldchanges","padj","include"]
    
    #print(df)
    df.index = df["names"]
    df = df.loc[df["include"].notna()]
    plt.figure(figsize=(5,10))
    ax=plt.subplot(1,1,1)

    #plt.scatter(df["logfoldchanges"],-1*np.log(df["padj"]))
    grey_rows = df.loc[(df["padj"] >pmin)].index
    plt.scatter(df.loc[grey_rows,"logfoldchanges"],-1*np.log10(df.loc[grey_rows,"padj"]),c="grey")

    middle_rows = df.loc[((df["logfoldchanges"] < lmax) | (df["logfoldchanges"] >lmin)) &(df["padj"] <pmin)].index
    plt.scatter(df.loc[middle_rows,"logfoldchanges"],-1*np.log10(df.loc[middle_rows,"padj"]),c="b")

    good_rows = df.loc[((df["logfoldchanges"] > lmax) | (df["logfoldchanges"] <lmin)) &(df["padj"] <pmin)].index
    plt.scatter(df.loc[good_rows,"logfoldchanges"],-1*np.log10(df.loc[good_rows,"padj"]),c="r")

    #print(len(good_rows))
    for gene in good_rows:   
        xy = (df.loc[gene,"logfoldchanges"], -1*np.log10(df.loc[gene,"padj"]))
        ax.annotate(gene, xy=xy, textcoords='data') # <--
    ylim = ax.get_ylim()
    xlim = ax.get_xlim()
    plt.plot([lmin,lmin],list(ax.get_ylim()),"--",c="k")
    plt.plot([lmax,lmax],list(ax.get_ylim()),"--",c="k")
    plt.plot(list(ax.get_xlim()),[-1*np.log10(pmin),-1*np.log10(pmin)],"--",c="r")
    #print(-1*np.log(pmin))
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    plt.xlabel("logfoldchanges")
    plt.ylabel("-logpadj")

def save_filtered_rankgenesgroups(adata_no_doublets,directory):
    names = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]["names"])
    logchnges = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]["logfoldchanges"])
    pvals_adj = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pvals_adj'])
    pvals = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pvals'])
    pts = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pts'])
    #pts_rest = pd.DataFrame(adata_no_doublets.uns["rank_genes_groups_filtered"]['pts_rest'])
    method = adata_no_doublets.uns["rank_genes_groups"]['params']['method']
    for i in names.columns:
        names_include = names[i].notna()
        #print(names_include)
        n_df = pd.DataFrame(index=names[i][names_include], columns=["logfoldchanges","pvals_adj","pvals","pt_expressing"])#,"pt_other_expressing"])
        n_df["logfoldchanges"] = logchnges[i].values[names_include]
        n_df["pvals_adj"] = pvals_adj[i].values[names_include]
        n_df["pvals"] = pvals[i].values[names_include]
        n_df["pt_expressing"]=pts.loc[names[i][names_include],i].values
        #n_df["pt_other_expressing"]=pts_rest.loc[names[i][names_include],i].values
        if "/" in i:
            i=i.replace("/","_")
        
        n_df.to_csv(directory+"filtered_"+method+"_all_cells_"+i+".csv")
def boxplot_sample_proportions(adata, x_value, color_value,hue="treatment",figsize=(10,5), plottype="box",order=None,hue_order=None,edgecolor=False,swap=False):
    tmp = adata.obs.groupby([x_value,color_value])[color_value].count().unstack(color_value).fillna(0)

    m=tmp.divide(tmp.sum(axis=1), axis=0)
    props = []

    i=0
    if hue+"_colors" in adata.uns and not swap:
        color_dict = dict(zip(adata.obs[hue].cat.categories,adata.uns[hue+"_colors"]))
    elif color_value+"_colors" in adata.uns and swap:
        color_dict = dict(zip(adata.obs[color_value].cat.categories,adata.uns[color_value+"_colors"]))
    else:
        color_dict=None
    for sample in m.index:

        for celltype in m.columns:
            vals = [sample,m.loc[sample,celltype],celltype,adata.obs.loc[adata.obs[x_value]==sample,hue].unique()[0]]
            props.append(vals)
            i+=1
    props_df = pd.DataFrame(props,columns=[x_value,x_value+"_proportion",color_value,hue])
    #sns.boxplot(x="celltype", y="sample_proportion", hue="treatment", data=tips)
    props_df[hue]=props_df[hue].astype("category")
    plt.figure(figsize=figsize)
    if swap:
        old_hue = hue
        hue=color_value
        color_value=old_hue
        old_hue_order=hue_order
        hue_order=order
        order=old_hue_order
    if plottype=="box":

        p=sns.boxplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df, palette=color_dict,hue_order=hue_order,linewidth=3)
        
        if edgecolor==True:
            for i,box in enumerate(p.artists):
                box.set_edgecolor(box.get_facecolor())
                r,g,b,a = box.get_facecolor()
                box.set_facecolor((r,g,b,.3))
            swarm_palette=color_dict
        else:
            swarm_palette={i:"white" for i in color_dict}
        if hue_order is None:
            hue_order = adata.obs[hue].cat.categories
        sns.swarmplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df,  dodge=True,palette=swarm_palette,edgecolor="white",linewidth=.7,size=3.2, hue_order=hue_order)
        plt.legend(p.artists,hue_order)
    if plottype=="bar":
        p=sns.barplot(x=color_value, y=x_value+"_proportion", hue=hue, data=props_df,palette=color_dict,hue_order=hue_order)
    p.set_xticklabels(p.get_xticklabels(),rotation=90)



def gsea_pairwise(adata, col, direction = "greater", minlfc = 1.0, maxpadj = .05, gene_sets = ["KEGG_2019_Human","MSigDB_Hallmark_2020","HMDB_Metabolites","WikiPathways_2019_Human","GO_Biological_Process_2018"],title=None):
    import gseapy as gp
    # options for gene sets https://maayanlab.cloud/Enrichr/#stats
    df=pd.concat([pd.DataFrame(adata.uns["rank_genes_groups"]["names"])[col],pd.DataFrame(adata.uns["rank_genes_groups"]["logfoldchanges"])[col],pd.DataFrame(adata.uns["rank_genes_groups"]["pvals_adj"])[col],pd.DataFrame(adata.uns["rank_genes_groups_filtered"]["names"])[col]],axis=1)
    df.columns=["names","logfoldchanges","padj","include"]
    df.index = df["names"]
    df = df.loc[df["include"].notna()]
    if direction == "greater":
        enr = gp.enrichr(gene_list=list(df.loc[(df["logfoldchanges"]>minlfc)&(df["padj"]<maxpadj)]["names"].values),
                         gene_sets=gene_sets, cutoff=0.1)
    else:
        enr = gp.enrichr(gene_list=list(df.loc[(df["logfoldchanges"]<minlfc)&(df["padj"]<maxpadj)]["names"].values),
                         gene_sets=gene_sets, cutoff=0.1)

    if title is None:
        title = col
    gp.dotplot(enr.results, title=title, figsize=(4,8),  top_term=20)
    return enr
