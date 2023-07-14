import scanpy as sc
import anndata
import os
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
import umap
#import scrublet as scr
import importlib
import sys
#sys.path.insert(0,'/data/cb/nyquist/breast_milk/breastMilk/')
import plotting_helpers as hh
sys.path.insert(0, '/data/cb/nyquist/cd8_deplete/')
import condition_plotting_helpers as cph
#import cd8_helpers as cd8
import seaborn as sns


LOCs = {"LOC101925857":"TRAC",
             "LOC102144039":"Mamu-DRA",
                  "LOC102122418":"CYBB*",
                        "LOC102136468":"MafaHLA-DPB1",
                              "LOC102140945":"FCGR3",
                                    "LOC102147203":"GBP1",
       "LOC102142617":"HLA1-B14A",
                "LOC102129434":"RPS11-like",
                "LOC107129205":"TRGV108B",
                "LOC102115168":"TRDC",
                "LOC102143603":"PWWP-pseudo",
                "LOC107130355":"TRAV-HPB-MLT",
                "LOC102115251":"KIR3DL2",
                "LOC107131025":"IGLV1-BL2",
                "LOC102144764":"KIR3DL3",
                "LOC102116611":"KIR2DL3",
                "LOC102139636":"HLA-DB1",
                "LOC102132533":"KLRC2",
                "LOC102140945":"FCGR3A",
                "LOC102132169":"KLRC2-like",
                "LOC102122509":"KLRC4",
                "LOC102115805":"IGHV-MPC11",
                "LOC102133674":"KLRC2-like.2",
                "LOC102116023":"MAP3K8",
                "LOC102134776":"CYB5A",
                "LOC102145001":"SCART1",
                "LOC102124348":"uncharacterized gene",
                "LOC102128672":"TRG2C",
                "LOC102133485":"IFITM3-pseudo",
                "LOC102142127":"uncharacterized gene.2",
                "LOC102145938":"IFITM3-like",
                "LOC102134129":"IFI27-pseudo",
                "LOC102123623":"HIRA",
                "LOC102139778":"uncharacterized gene.3",
                "LOC102136964":"LAIR1",
                "LOC102135811":"BOLA2-pseudo",
                "LOC102143485":"uncharacterized gene.4",
                "LOC102116308":"IGHM-like",
                "LOC102132314":"H2AX",
                "LOC102142866":"FANCD2",
                "LOC102141578":"H3C1-like",
                "LOC102119226":"H1-5",
                "LOC107128728":"uncharacterized gene.5",
                "LOC107126576":"uncharacterized gene.6",
                "LOC102135858":"HGMB2-pseudo",
                "LOC102121131":"H1-3",
                "LOC102123264":"H1-4",
                "LOC102128802":"TUBA1C",
                "LOC102138750":"H1-2",
                "LOC102137365":"H2AC12",
                "LOC102136402":"APITD1",
                "LOC102127357":"APOBEC3D",
                "LOC102146491":"H2AZP1",
                "LOC102134721":"ANP32E-pseudo",
                "LOC107129707":"uncharacterized gene.7",
                "LOC102118091":"H2AC14",
                "LOC102124147":"SPIB",
                "LOC102131515":"LRR1",
                "LOC102144925":"uncharacterized gene.8",
                "LOC102134772":"DEK-pseudo",
                "LOC102142158":"AKR7A2",
                "LOC102119092":"H4",
                "LOC102133926":"HGMB1-like",
                "LOC102134966":"COA4-h",
                "LOC102139881":"MAP3K20",
                "LOC102129658":"uncharacterized gene.9",
                "LOC102128706":"IGHA2-like",
                "LOC102142558":"IGHG1-like",
                "LOC102122028":"RPL40-pseudo",
                "LOC102120885":"RPLP0-pseudo",
                "LOC102121775":"Gogo-B*0103A-like",
                "LOC102131664":"HLA-A-24-like",
                "LOC102137590":"RPLP0-pseudo.2",
                "LOC102136805":"CCL4",
                "LOC102122037":"GOLIM4-like",
                "LOC102142222":"HLA-B-15-like",
                "LOC101867070":"uncharacterized gene.10",
                "LOC102136862":"HLA-DPA1-like",
                "LOC102139771":"RPL14-pseudo",
                "LOC102123006":"ARL6IP1-pseudo",
                "LOC102141835":"HLA-B7-like",
                "LOC102116527":"Gogo-B*0103A-like.2",
                "LOC102140002":"Gogo-B*0103A-like.3",
                "LOC102132418":"EEF1A1-pseudo",
                "LOC102143352":"RPL40-S18-pseudo",
                "LOC102137002":"uncharactarized gene.11",
                "LOC101865988":"uncharactarized gene.12",
                "LOC107128653":"uncharactarized gene.13",
                "LOC102122467":"uncharactarized gene.14",
                "LOC107129085":"uncharactarized gene.15",
                "LOC101925240":"HSPA1A",
                "LOC102119366":"CCL4-like",
                "LOC102125474":"PKM2-pseudo",
                "LOC102116897":"MafaHLA-B7",
                "LOC102128579":"ARIH2-like",
                "LOC102116151":"HLA-B-37-like",
                "LOC102141029":"RPL21",
                "LOC107126798":"uncharactarized gene.16",
                "LOC102143405":"Mafa-B-PATR-B-like",
                "LOC102142031":"ZNF33A-like",
                "LOC102131261":"BIN1",
                "LOC107130579":"uncharactarized gene.17",
                "LOC102133096":"SET-like",
                "LOC102139567":"LRRC37A2",
                "LOC102142617":"HLA1-B14A",
                "LOC102129434":"40S S11",
                "LOC107129205":"TRGV108B",
                "LOC102136846":"HBA1/2/3",
                           "LOC102136192":"HBA1/2/3-like",
                           "LOC102145413":"Mafa-DOB",
                           "LOC102134487":"CYB561A3",
                           "LOC102133296":"IGKV1-39-like",
                           "LOC102140918":"SH3BP5",
                           "LOC102119202":"PPP1R16B",
                           "LOC102126451":"TNFRSF13B",
                           "LOC102133073":"AKAP2",
                           "LOC107130289":"IGHG3-like",
                           "LOC102144414":"Mafa-DRB1-13-like",
                           "LOC107128164":"uncharactarized gene.18",
                           "LOC102126820":"GRAP",
                           "LOC102141176":"Mafa-DQA1-like",
                           "LOC102144867":"RPS11-like",
                           "LOC107128718":"uncharactarized gene.19",
                           "LOC102122749":"NAPSA-like",
                           "LOC102144782":"Mafa-DQB1-like",
                           "LOC102129412":"LIMD1",
                           "LOC102137481":"Mafa-DOA",
                           "LOC107126411":"IGLC6-like",
                           "LOC102144295":"HLA-A-11-like",
                           "LOC102127115":"RPS23-pseudo",
                           "LOC102133073":"AKAP2",
                           "LOC102136305":"uncharactarized gene.20",
                           "LOC102146935":"IFITM3-like",
                           "LOC102130734":"A2M",
                           "LOC102142030":"TM4SF1-like",
                           "LOC102119459":"PPFIBP1",
                           "LOC102140918":"SH2BP5",
                           "LOC102146822":"PTPRJ", "LOC102124330":"CXCL5","LOC102120556":"STAT5A",
                           "LOC102130733":"IFNL3","LOC102138968":"IFNL2","LOC102147203":"GBP1","LOC102146847":"GBP3",
                            "LOC102134512":"CCL23","LOC102134112":"CCL15","LOC102125111":"PF4","LOC102116354":"ZNF675",
                           }


def run_densmap(adata):
    '''
    Note: this requires that you change the scanpy code, so if you update scanpy this will stop working 
    '''
    copy_old_umap = False
    if "X_umap" in adata.obsm:
        umap_tmp = adata.obsm["X_umap"]
        copy_old_umap = True
    sc.tl.umap(adata,densmap=True,densmap_kwds={"graph_dists":adata.obsp['distances'],"lambda":2.0,"var_shift":0.1,"frac":0.3})
    adata.obsm["X_densmap"] = adata.obsm["X_umap"]
    if copy_old_umap:
        adata.obsm["X_umap"] = umap_tmp

def clustering_qc_plots(gran_adata,clusters_ordered_by_celltype,cluster_name):
    import matplotlib.gridspec as gridspec
    from matplotlib import ticker

    plt.rcParams.update({'font.size': 15})
    # Plot figure with subplots of different sizes
    fig = plt.figure(figsize=(35,15))

    # set up subplot grid
    grid=gridspec.GridSpec(1,4, width_ratios=(.8,.8,1.2,1))
    gs1 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = grid[0]) # clusters and assignments
    gs2 = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec = grid[1]) # monkeys and samples
    gs3 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = grid[2], wspace=1) # marker genes
    gs4 = gridspec.GridSpecFromSubplotSpec(6, 1, subplot_spec = grid[3]) # composition and QC

    ax_leiden = fig.add_subplot(gs1[0, 0])
    sc.pl.umap(gran_adata, color="leiden",show=False,ax=ax_leiden,legend_loc="on data")
    ax_celltypes = fig.add_subplot(gs1[1, 0])
    sc.pl.umap(gran_adata, color=cluster_name,show=False,ax=ax_celltypes,legend_loc="on data")

    ax_monkey = fig.add_subplot(gs2[0, 0])
    sc.pl.umap(gran_adata, color="M.Number",show=False,ax=ax_monkey)
    ax_monkey.get_legend().remove()
    ax_monkey.set_title("Monkey Number")
    ax_sample = fig.add_subplot(gs2[1, 0])
    sc.pl.umap(gran_adata, color="sample",show=False,ax=ax_sample)
    ax_sample.get_legend().remove()

    ax_mg = fig.add_subplot(gs3[0, 0])
    ax_out=sc.pl.rank_genes_groups_dotplot(gran_adata,groupby="leiden", categories_order=clusters_ordered_by_celltype,groups=clusters_ordered_by_celltype,
                                    n_genes=3,dendrogram=False,show=False,swap_axes=True,ax=ax_mg)
    replace_rankgenesgroups_dotplot_locs(ax_out)
    qcplots(gran_adata,groupby="leiden",gs4=gs4,fig=fig)
    return gs1,gs2,gs2,gs4

def qcplots(gran_adata, groupby="leiden", gs4=None,fig=None, donor_colname = "M.Number",sample_colname="sample",include_stackedbars=True):
    import matplotlib.gridspec as gridspec
    from matplotlib import ticker
    if gs4 is None:
        if include_stackedbars:
            fig=plt.figure(figsize=(7,15))
            gs4 = gridspec.GridSpec(6,1)
        else:
            fig=plt.figure(figsize=(7,11))
            gs4 = gridspec.GridSpec(4,1)
    ax_tc = fig.add_subplot(gs4[0, 0])
    #else:
        #gs4 = ax.get_subplotspec()
        #ax_tc=ax
    sc.pl.violin(gran_adata, "total_counts",groupby=groupby,rotation=90,ax=ax_tc,show=False,stripplot=False)
    ax_tc.set_xlabel("")
    ax_tc.set_xticklabels([])
    ax_tc.set_xticks([])
    ax_tc.set_ylabel("n_UMI")
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1))
    ax_tc.yaxis.set_major_formatter(formatter)
    ax_mito = fig.add_subplot(gs4[1, 0])
    sc.pl.violin(gran_adata, "percent_mito",groupby=groupby,rotation=90,ax=ax_mito,show=False, stripplot=False)
    ax_mito.set_xlabel("")
    ax_mito.set_xticklabels([])
    ax_mito.set_xticks([])
    ax_mito.set_ylabel("%mito")
    ax_genes = fig.add_subplot(gs4[2, 0])
    sc.pl.violin(gran_adata, "n_genes_by_counts",groupby=groupby,rotation=90,ax=ax_genes,show=False, stripplot=False)
    ax_genes.set_xlabel("")
    ax_genes.set_xticklabels([])
    ax_genes.set_xticks([])
    ax_genes.set_ylabel("n_genes")
    formatter_g = ticker.ScalarFormatter(useMathText=True)
    formatter_g.set_scientific(True) 
    formatter_g.set_powerlimits((-1,1))
    ax_genes.yaxis.set_major_formatter(formatter_g)
    ax_doublet = fig.add_subplot(gs4[3, 0])
    sc.pl.violin(gran_adata, "doublet_scores",groupby=groupby,rotation=90,ax=ax_doublet,show=False, stripplot=False)
    ax_doublet.set_ylabel("doublet\nscores")
    if include_stackedbars:
        ax_doublet.set_xlabel("")
        ax_doublet.set_xticklabels([])
        ax_doublet.set_xticks([])
        ax_sample = fig.add_subplot(gs4[4, 0])
        hh.normalized_stacked_bar_plot(gran_adata, groupby,sample_colname,ax=ax_sample,legend=False)
        ax_sample.set_xlabel("")
        ax_sample.set_xticklabels([])
        ax_sample.set_xticks([])
        ax_sample.set_ylabel("pct cells")
        ax_monkey = fig.add_subplot(gs4[5, 0])
        hh.normalized_stacked_bar_plot(gran_adata, groupby,donor_colname,ax=ax_monkey,legend=False)
        ax_monkey.set_ylabel("pct cells")

def replace_rankgenesgroups_dotplot_locs(ax_out):
    lbls = []
    for label in ax_out["mainplot_ax"].get_yticklabels():
        if label.get_text() in LOCs:
            label.set_text(LOCs[label.get_text()]+"*")
        lbls.append(label)
    ax_out["mainplot_ax"].set_yticklabels(lbls,fontdict={'fontsize':10})
    lbls = []
    for label in ax_out["mainplot_ax"].get_xticklabels():
        if label.get_text() in LOCs:
            label.set_text(LOCs[label.get_text()]+"*")
        lbls.append(label)
    ax_out["mainplot_ax"].set_xticklabels(lbls,fontdict={'fontsize':10})

    
def save_pseudobulk_counts(count_adata,processed_metadata=None,celltype_col=None,sample_col="sample"):
    '''
    count_adata should be an adata object with raw counts saved in the 'X' slot
    processed_metadata [a dataframe] is probably the obs of whatever adata object you have been working with, but if it is None, use the metadata from count_adata
    celltype_col is the column in the processed_metadata dataframe on which you would like to split your sample pseudobulks. If it is None, a single pseudobulk is created for each sample.
    sample_col is the other column to split your cells on (usually just called sample...)
    '''
    if processed_metadata is None:
        processed_metadata = count_adata.obs

    exp = pd.DataFrame(count_adata.X,index=count_adata.obs_names,columns=count_adata.var_names).loc[processed_metadata.index]
    exp[sample_col] = processed_metadata[sample_col]

    if celltype_col is not None:
        exp[celltype_col] = processed_metadata[celltype_col]
        pseudobulk = exp.groupby([celltype_col,sample_col]).sum()
        pseudobulk[str(celltype_col)+"_"+str(sample_col)] = [i[0]+"_"+i[1] for i in pseudobulk.index]
        pseudobulk.index = pseudobulk[str(celltype_col)+"_"+str(sample_col)]
        pseudobulk.columns = map_loc_list(pseudobulk.columns)
        return pseudobulk[[c for c in pseudobulk.columns if c not in [sample_col,celltype_col,str(celltype_col)+"_"+str(sample_col)]]]
    else:
        return exp.groupby([sample_col]).sum()

def save_metadata_for_pseudobulk_deseq(processed_metadata, celltype_col = None, sample_col="sample",cell_counts=True):
    meta_cols = processed_metadata.columns
    if celltype_col is not None:
        grp_name = "TMP sample cols"
        processed_metadata[grp_name] = [processed_metadata.loc[i,celltype_col]+"_"+processed_metadata.loc[i,sample_col] for i in processed_metadata.index]
        df_index = list(processed_metadata[grp_name].unique())
    else:
        df_index = list(processed_metadata[sample_col].unique())
        grp_name = sample_col
    meta = pd.DataFrame(index=df_index, columns=meta_cols)

    for c in meta_cols:
        processed_metadata[c] = processed_metadata[c].astype(str)

        v = processed_metadata.groupby(grp_name)[c].unique()
        meta[c] = meta.index.map({i:v.loc[i][0] for i in v.index})
    if cell_counts:
        cnts = processed_metadata.groupby([grp_name]).count().iloc[:,0]
        meta["Cell_Number"]=cnts.loc[meta.index]
    return meta


def condition_split_violin(adata,cat,vertical_split,genenames,vertical_order = None, violin_order=None):
    if isinstance(genenames,dict):
        genes=list(genenames.keys())
    elif isinstance(genenames,list):
        genes = genenames
        genenames = dict(zip(genes,genes))
    else:
        genes = [genenames]
        genenames = {genes[0]:genes[0]}
    if vertical_order is None:
        vertical_order = adata.obs[vertical_split].cat.categories
    n_rows = len(adata.obs[vertical_split].unique())
    fig,ax = plt.subplots(n_rows,len(genes),figsize=(len(genes)*len(adata.obs[cat].unique())/2.0,n_rows*1.5),sharey=True,gridspec_kw={"hspace":0,"wspace":0})
    if cat+"_colors" in adata.uns:
        palette = dict(zip(adata.obs[cat].cat.categories,adata.uns[cat+"_colors"]))
    else:
        palette = "Blues"
    if violin_order is None:
        violin_order = adata.obs[cat].cat.categories
    for g_i,gene in enumerate(genes):
        for i,t in enumerate(vertical_order):
            igg= adata[adata.obs[vertical_split]==t]

            df=pd.DataFrame(igg.raw[:,gene].X,index=igg.obs_names, columns=[gene,])
            df[cat]= igg.obs[cat]
            sns.violinplot(y=gene,x=cat,data=df,scale="width",ax=ax[i][g_i],palette=palette,order=violin_order)
            if i <= len(ax)-2:
                ax[i][g_i].set_xticklabels([])
                ax[i][g_i].set_xticks([])
            ax[i][g_i].set_xlabel("")
            if g_i > 0:
                ax[i][g_i].tick_params(left=False,labelleft=False)
                #plt.setp(ax[i][g_i].get_yticks(), visible=False)
                ax[i][g_i].set_ylabel("")
            else:
                ax[i][g_i].set_ylabel(t,rotation=0,labelpad=5+len(t)*2,size=10)
        ax[0][g_i].set_title(genenames[gene])
        for tick in ax[i][g_i].get_xticklabels():

            tick.set_rotation(90)



def colors_from_spreadsheet(sheet_name,celltype_col="val",value_col='color'):
    import gspread
    from oauth2client.service_account import ServiceAccountCredentials
    SCOPES = ['https://www.googleapis.com/auth/spreadsheets','https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name('/data/cb/scratch/nyquist/deeper_sequencing_analysis/token.json', SCOPES)
    client = gspread.authorize(creds)
    sheet = client.open('CD8 Depletion Paper Color Schemes and Celltype Names')
    df = pd.DataFrame.from_dict(sheet.worksheet(sheet_name).get_all_records())
    df[celltype_col]=df[celltype_col].astype(str)
    df.index=df[celltype_col]

    return df[value_col].to_dict()

def correct_all_adata_colors(adata, general_celltype_names = ["General Celltypes"],major_nkt_names = ["Major NK/T Clusters"], nkt_subclusters_name=["NK/T Subclusters"],treatment_name=["treatment"],mnumber_name=["M.Number"]):
    columns_to_obs = {"General Celltypes":general_celltype_names,
                        "Major NK/T Clusters":major_nkt_names,
                        "NK/T Subclusters":nkt_subclusters_name,
                        "treatment":treatment_name,
                        "M.Number":mnumber_name}
    for sheet_col,adata_cols in columns_to_obs.items():

        col_dict = colors_from_spreadsheet(sheet_col)
        for col in adata_cols:
        
            if col in adata.obs:
                hh.set_colors_from_dict(adata,col_dict,col)
def update_names(adata, nkt_sub_col = ["NK/T Subclusters"], old_col="val",new_col="new_name"):
    columns_to_obs = {"NK/T Subclusters": nkt_sub_col}
    for sheet_col,adata_cols in columns_to_obs.items():
        col_dict = colors_from_spreadsheet(sheet_col, celltype_col=new_col, value_col="color")
        name_dict = colors_from_spreadsheet(sheet_col, celltype_col=old_col, value_col=new_col)
        for col in adata_cols:
            if col in adata.obs:
                adata.obs[col] = adata.obs[col].astype(str)
                adata.obs[col]=adata.obs[col].map(name_dict)
                hh.set_colors_from_dict(adata, col_dict, col)

def save_filtered_rankgenesgroups_googledrive(adata_no_doublets,sheetname):
    import gspread
    from oauth2client.service_account import ServiceAccountCredentials
    SCOPES = ['https://www.googleapis.com/auth/spreadsheets','https://www.googleapis.com/auth/drive']
    creds = ServiceAccountCredentials.from_json_keyfile_name('token.json', SCOPES)
    client = gspread.authorize(creds)
    try:
        sheet = client.open(sheetname)
    except gspread.SpreadsheetNotFound as e:
        sheet = client.create(sheetname)
        print("created new spreadsheet")
    client.insert_permission(sheet.id,"sarahnyquist@gmail.com",perm_type="user",role="owner",notify=False)
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
        n_df.index = map_df_index_locs(n_df.index)
        if i in [s.title for s in sheet.worksheets()]:
            w = sheet.worksheet(i)
            w.clear()
        else:
            w = sheet.add_worksheet(title=i, rows=str(n_df.shape[0]+1),cols=str(n_df.shape[1]+1))
        df_to_worksheet(n_df, w)
def map_loc_list(genes):
        loc_mapped = []
        for n in genes:
            if n in LOCs:
                loc_mapped.append(LOCs[n]+"*")
            else:
                loc_mapped.append(n)
        return loc_mapped

def df_to_worksheet(df,worksheet):
    worksheet.resize(df.shape[0]+1,df.shape[1]+1)
    old_cols = df.columns
    df[df.index.name] = df.index
    col_order = [df.index.name]+list(old_cols)
    l = df[col_order].values.tolist()
    lol = [col_order]+l
    worksheet.update('A1:E'+str(len(lol)),lol)

