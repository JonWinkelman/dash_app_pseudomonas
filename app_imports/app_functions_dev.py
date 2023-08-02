from Bio import Phylo
import pandas as pd
import os
import json
import plotly.graph_objs as go
import numpy as np
from jw_utils import jw_draw_tree as jdt
# from ete3 import ncbi_taxonomy
# ncbi_tax = ncbi_taxonomy.NCBITaxa()
from jw_utils import ncbi_datasets_fxs as nfx
from jw_utils import jw_ncbi_taxonomy as jnt

def acc2name(acc2name_path):
    """Return accession_to_name dict from AssemblyAccession_to_SpeciesName.json in summary data""" 
    with open(acc2name_path, 'r') as f:
        return json.load(f)
    
def name2acc(acc2name_path):
    """Return name to accession dict"""
    d = acc2name(acc2name_path)
    return {key:val for val, key in d.items()}


def gene_protein_dict(path_to_gff3, gene_name):
    'get all features that have user-input gene as their parent'
    with open(path_to_gff3, 'r') as f:
        gene_prot_dict = {}
        att_dict = {}
        for line in f:
            if line.find(gene_name) != -1 and line.find('Parent=gene') != -1:
                attribute_lst= line.split('\t')[-1].split(';')
                for ele in attribute_lst:
                    key_val = ele.split('=')
                    att_dict[key_val[0]] = key_val[1]
                if not gene_prot_dict.get(gene_name):
                    gene_prot_dict[gene_name] = {'protein':[att_dict['ID']],
                                                 'product':[att_dict['product']]}
                else:
                    gene_prot_dict[gene_name]['protein'].append(att_dict['ID'])
                    gene_prot_dict[gene_name]['product'].append(att_dict['product'])
    return gene_prot_dict



def make_arrow_trace(protein, par, strand, prod, start, end, br, HOG, **kwargs): #color=None,
    if kwargs['color']:
        color = kwargs['color']
    else:
        color = 'rgb(100,100,100)'

    'return a trace with an arrow drawn for the input feature'
    
    trace = go.Scatter(
            x = [start, start, br, br, end,br,br,start,start],
            y = [0, 0.25,0.25,0.5,0,-0.5,-0.25,-0.25,0],
            mode = 'lines',
            fill='toself',
            line = {'width': 1.5,
                    'color':color},
            name = protein,
            text = f'{protein}<br><br>{par}<br>strand: {strand}<br>{prod}<br>start: {start}<br>end: {end}<br>HOG: {HOG}',
            hoverinfo = 'text')
    return trace



def make_map(trimmed_df, feature_name = None, height = 150, yrange = [-1,1], colors=None, bgcolor='rgb(250,250,250)'):
    traces = []
    xrange = [0,0]
    #make traces for each feature
    for i,protein in enumerate(trimmed_df.index):
        #hoverinfo variables
        par =  trimmed_df.loc[protein,'Parents']
        prod = trimmed_df.loc[protein,'products']
        start = trimmed_df.loc[protein, 'starts']
        end = trimmed_df.loc[protein, 'ends']
        strand = trimmed_df.loc[protein,'strand']
        #orthogroup = trimmed_df.loc[protein,'Orthogroups']
        HOG = trimmed_df.loc[protein,'HOGs']
        #make backbone trace between features
        if i< (len(trimmed_df.index)-1):
            next_orf_start = trimmed_df.iloc[(i+1),0]
            traces.append(go.Scatter(x = [end, next_orf_start],
                          y = [0,0],
                          mode = 'lines',
                          line = {'width': 3,
                                  'color':'black'},
                          showlegend = False,
                          hoverinfo = None))
        if strand == '-':
            start =  end
            end = trimmed_df.loc[protein, 'starts']   
        l = (end-start)     #lenght of the arrow 
        br = start+(0.65*l) #defines where head of arrow is placed
         #make feature traces, highlighting the feature of interest
        if protein == feature_name:
            traces.append(make_arrow_trace(protein, par, strand, prod, start, end, br, HOG, color=colors['t_red']))
            xrange  =[start-10000, end+10000]
        else:
            traces.append(make_arrow_trace(protein, par, strand, prod, start, end, br, HOG, color='grey'))
        
    dl = {'data':traces,
            'layout':go.Layout(#title = f'local genome map around {feature_name}',
                               paper_bgcolor='rgb(255,255,255)',
                               plot_bgcolor=bgcolor,
                               width = 700,
                               height = height,
                               margin={'t':0, 'b':0, 'l':0, 'r':0},
                               showlegend = False,
            )} 
    fig = go.Figure(dl)  
    fig.update_yaxes(showgrid=False, showticklabels=False)
    fig.update_yaxes(range=yrange)
    fig.update_xaxes(range=xrange)
    return fig
    





def create_HOGgene_tree(tree, selected_name, ingroup, outgroup, orthoparser_obj=None, colors=None):

    x_coords = jdt.get_x_coordinates(tree)
    y_coords = jdt.get_y_coordinates(tree)
    line_shapes = []
    jdt.draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
    ##replace_periods in the accesions --- orthofinder and other prgrms als do this
    accession_list_edit = [acc.replace('.','_') for acc in list(orthoparser_obj.accession_to_name.keys())]
    name_gene_dict = {}
    node_size_dict = {}
    node_color_dict = {}
    for i,cl in enumerate(my_tree_clades):
        if cl.name:
            node_size_dict[cl.name] = 4
            node_color_dict[cl.name] = '#92A0A9' #grey'
            name_gene_dict[cl.name] = cl.name
            for accession, accession_ed in zip(list(orthoparser_obj.accession_to_name.keys()),accession_list_edit):
                if cl.name.startswith(accession_ed):
                    ind = len(accession_ed)
                    name_gene_dict[cl.name]  = orthoparser_obj.accession_to_name[accession] + '  ' + cl.name[ind+1:]
        else:
            node_size_dict[str(i)] = 2
            node_color_dict[str(i)] = '#92A0A9' #grey'
            name_gene_dict[str(i)] = 'N'+ str(i) 
            
    for accession in ingroup:
        accession = accession.replace('.','_')
        for key in name_gene_dict.keys():
            if key.startswith(accession):
                node_color_dict[key] = colors['t_green']  #seagreen trestle
                node_size_dict[key] = 10

                
    for accession in outgroup:
        accession = accession.replace('.','_')
        for key in name_gene_dict.keys():
            if key.startswith(accession):
                node_color_dict[key] = colors['t_blue'] #?? orange
                node_size_dict[key] = 10
                
    for key in name_gene_dict.keys():
        if key:
            selected_name_edit = selected_name.replace('.', '_')
        if key.startswith(selected_name_edit):
            node_size_dict[key] = 12
            node_color_dict[key] = colors['t_red']
              
    axis = dict(showline=False,
                zeroline=False,
                showgrid=False,
                showticklabels=False,
                title=''  # y title
        )

    data = dict(type='scatter',
                x=X,
                y=Y,
                mode='markers',
                marker=dict(color=list(node_color_dict.values()),
                            size =list(node_size_dict.values())),
                text=list(name_gene_dict.values()),  # vignet information of each node
                hoverinfo='text',
        )
  
    layout = dict(#title='Gene tree: Genes in orthogroup '+ selected_orthogroup,
                  paper_bgcolor='rgba(0,0,0,0)',
                  #dragmode="select",
                  font=dict(family='Balto', size=14),
                  #width=750,
                  height=500,#len(tree.get_terminals())*8,
                  autosize=True,
                  showlegend=False,
                  xaxis=dict(showline=False,
                             zeroline=False,
                             showgrid=False,  # To visualize the vertical lines
                             ticklen=4,
                             showticklabels=False,
                             title=''),
                  yaxis=axis,
                  hovermode='closest',
                  shapes=line_shapes,
                  plot_bgcolor='rgb(250,250,250)',
                  legend={'x': 0, 'y': 1},
                  margin={'b': 0, 'l': 20, 'r': 0, 't': 0}
        )
    fig = dict(data=[data], layout=layout)   
    return fig


def gene_to_orthogroup(path_to_orthogroups, species_accession):
    'return a dict (genes:orthogroups) mapping genes in a gven species to its orthogroups'
    df = pd.read_csv(path_to_orthogroups, sep = '\t',low_memory=False)
    df = df.loc[:,[species_accession, 'Orthogroup']]
    orth_dict = {}
    for ind in df.index:
        if type(df.iloc[ind, 0]) == str:
            orth_dict[df.iloc[ind, 1]] = df.iloc[ind, 0].split(',')
    #make new dict
    gene_to_orthogroup_dict = {}
    for orth in orth_dict.keys():
        for val in orth_dict[orth]:
            gene_to_orthogroup_dict[val] = orth
    return gene_to_orthogroup_dict



def get_df_for_feature(dop_obj, assembly_accession, feature, fts=15):
    """Return a local df containing features near the passed feature
    
    fts (int): num features to show in map on each side of of feature of interest 
    """

    df = dop_obj.make_genome_annotation_df(assembly_accession)
    if feature in df.index:       
        f_index = df.index.get_loc(feature)
        if f_index>=fts and f_index<= (len(df.index)-fts):
            df = df.iloc[(f_index-fts):(f_index+fts),:]
        
        elif f_index-fts<0:
            df = df.iloc[:(f_index+fts),:]
        elif f_index + fts > df.shape[0]:
            df = df.iloc[(f_index-fts):,:]
    else:
        print(f'{feature} not in in the dataframe index, e.g. index = {df.index[0]}')
        df = df.iloc[:10,:]
    return df
    

def get_accession_from_genetree_click(geneTree_clickData, dop_obj):
    name = ' '.join(geneTree_clickData['points'][0]['text'].split(' ')[:-1]).strip() 
    return dop_obj.name_to_accession[name]


def get_protein_from_genetree_click(geneTree_clickData):
    return geneTree_clickData['points'][0]['text'].split(' ')[-1]


def accession_taxid_d(path_to_summary):
    """Return simplified dict {acc:taxid} from ncbi summary.json."""
    
    with open(path_to_summary, 'r') as f:
        summary_d = nfx.make_summary_dict(json.load(f))
    acc_tax_d ={}
    for acc in summary_d:
        acc_tax_d[acc] = summary_d[acc]['org']['tax_id']
    return acc_tax_d



def get_lineage_rank_dict(taxID, ncbi_tax):
    """Return a lineage dict {rank:taxid} from a single taxid
    
    ranks = 'no rank', 'superkingdom','family', 'genus', 'phylum', 'class', 'order', 'species', 
    
    taxid (int): ncbi taxonomic ID that can be parsed by the ete3 library
    """
    d = ncbi_tax.get_rank(ncbi_tax.get_lineage(taxID))
    return  {key:val for val, key in d.items()} 



def get_rank_sciname_for_each_leaf(tree, rank, path_to_summary, ncbi_tax):
    """Return a dict with a sci_name at a given rank for each leaf {leaf:sci_name}
    e.g. given rank phylum: return {GCF_######:'Proteobacteria'}
    """

    acc_tax_d =  accession_taxid_d(path_to_summary,)   
    leaf_rank_sciname = {}
    for leaf in tree.get_terminals():
        leaf_taxid = acc_tax_d[leaf.name]
        leaf_rank_taxid = get_lineage_rank_dict(leaf_taxid, ncbi_tax)
        leaf_rank_taxid = leaf_rank_taxid.get(rank)
        if leaf_rank_taxid:
            leaf_rank_sciname[leaf.name] = list(ncbi_tax.get_taxid_translator([leaf_rank_taxid]).values())[0]
        else:
            leaf_rank_sciname[leaf.name] = None
    return leaf_rank_sciname



def get_list_of_all_nodes(subtree, clades=None):
    """Get a list of all clade objects, including leaves, in a Bio.Phylo tree"""
    if not isinstance(subtree, (Phylo.Newick.Tree, Phylo.Newick.Clade)):
        raise TypeError(f'object entered needs to be of type {Phylo.Newick.Tree} or {Phylo.Newick.Clade}, you entered {type(subtree)}' )
        
    if clades is None:
        clades=[]
    for cl in subtree.root:
        clades.append(cl)
        if cl.is_terminal():
             clades.append(cl)
        else:
            get_list_of_all_nodes(cl, clades)
    return clades



def get_nodes_assoc_with_ranks(tree, rank, path_to_summary, ncbi_tax):
    """"""
    leaf_rank_scinames_d = get_rank_sciname_for_each_leaf(tree, rank, path_to_summary, ncbi_tax)
    node_dict = {}
    for cl in get_list_of_all_nodes(tree):
        #check to see if all leaves in clade have same sci_name, if so then clade is that sci_name
        leaves =  cl.get_terminals()
        sci_names = []

        for leaf in leaves:
            sci_names.append(leaf_rank_scinames_d[leaf.name])
        sci_names = list(set(sci_names))
        if len(sci_names) == 1:
            node_dict[cl.name] = sci_names[0]
    unique_scinames = set(list(node_dict.values()))
    rank_nodes = {n:[] for n in unique_scinames}
    for node, sci_name in node_dict.items():
        rank_nodes[sci_name].append(node)
    return rank_nodes


def get_anc_node_each_rank(tree, rank, path_to_summary, ncbi_tax):
    """"""
    rank_nodes = get_nodes_assoc_with_ranks(tree, rank, path_to_summary,ncbi_tax)
    mrca_clades = {}
    for rank_name, clades in rank_nodes.items():
        leaves = []
        for cl in clades:
            clade = tree.find_any(cl) 
            if clade.is_terminal():
                leaves.append(clade)
        mrca_clades[rank_name] = tree.is_monophyletic(leaves)
    return mrca_clades