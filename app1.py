#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 08:38:52 2022

@author: jonwinkelman
"""
import dash
from dash import dcc, html
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import pandas as pd
import numpy as np
import os
import json
from jw_utils import parse_gff as pgf
from app_imports import app_elements_dev as ae
from app_imports import app_functions_dev as afn
from jw_utils import jw_draw_tree as jdt
from jw_utils import jw_ncbi_taxonomy as jnt
from ete3 import ncbi_taxonomy
ncbi_tax = ncbi_taxonomy.NCBITaxa()
from orthofinder_utils import dash_ortho_parser as dop

# =============================================================================
# set lab_name variable
# =============================================================================
lab_name = 'Mukherjee Lab'
USERNAME_PASSWORD_PAIRS = [['Trestle', 'MukherjeeLab']]
# =============================================================================
# |^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|
# =============================================================================

path_to_species_tree = './data/Species_Tree/SpeciesTree_rooted.txt'
#path_to_species_tree_intNodeNames = './data/Species_Tree/SpeciesTree_rooted_node_labels.txt'
path_to_species_tree = './data/Species_Tree/SpeciesTree_rooted_node_labels.txt'
path_to_gene_trees = './data/Resolved_Gene_Trees'
path_to_orthogroups =  './data/Orthogroups/Orthogroups.tsv'
acc2name_path = './data/summary_data/AssemblyAccession_to_SpeciesName.json'
path_to_summary = './data/summary_data/summaries.json'
colors = {
    't_blue': 'rgba(0,102,153,255)',
    't_green': 'rgba(61,174,43,255)',
    't_red': 'rgb(255,20,20)',
    'seagreen':'#2c8d42',
    'orange':'#F9B257',
    'purple':'rgba(52, 30, 170, 0.25)',
    'lavender':'rgba(170, 30, 125, 0.25)',
    'greyish': 'rgba(80, 80, 69, 0.25)'}


app = dash.Dash(__name__)
#auth =dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)
server = app.server

#################################################################################################################################################
##############################################Functions##########################################################################################


def get_annot_df_for_orthology(dop_obj, accession):
    df = dop_obj.make_genome_annotation_df(accession)
    df.index = [ind[4:] for ind in df.index]
    df = df.reset_index()
    df = df.loc[df['HOGs'] != 'no HOG',:]
    df = df.drop('Orthogroups',axis = 1)
    df.columns = ['Protein', 'Start', 'End', 'Strand', 'Product', 'Gene', 'HOG']
    return df


def get_orthologs(orth_species_1, orth_species_2):
    tree = jdt.read_treefile(path_to_species_tree)
    mrca_clade = tree.common_ancestor(orth_species_1, orth_species_2)
    dop_obj= dop.DashOrthoParser('./data', mrca_clade.name)
    df_bait = get_annot_df_for_orthology(dop_obj, orth_species_1 )
    df_prey = get_annot_df_for_orthology(dop_obj, orth_species_2 )
    return df_bait.merge(df_prey, on='HOG').set_index('Protein_x')


def get_orthology_counts(df_merged):
    """"""
    x_counts = {protein:0 for protein in df_merged.index}
    y_counts = {protein:0 for protein in df_merged.loc[:,'Protein_y']}
    for protein in df_merged.index:
        if protein in df_merged.index:
            x_counts[protein] += 1
    for protein in df_merged['Protein_y']:     
        if protein in list(df_merged['Protein_y']):
            y_counts[protein] += 1
    orthology_type = []
    for protein_x, protein_y in zip(df_merged.index, df_merged['Protein_y']):
        ratio  = f'{x_counts[protein_x]} to {y_counts[protein_y]}'
        orthology_type.append(ratio)
    df_merged['Orth_counts_x:y'] = orthology_type
    return df_merged






def make_dopObj_from_click(clickData):
    """Get node info from clickData and return dop_obj from appropriate level of HOG."""
    tree =  jdt.read_treefile(path_to_species_tree)
    if not clickData:
        node_name ='N0'
    else:
        node_name = clickData['points'][0]['text']
        node_name = node_name.split('<br>')[0].split(' ')[1]
        int_node_names = [cl.name for cl in tree.get_nonterminals()]
        if node_name not in int_node_names:
            node_name ='N0'
    return dop.DashOrthoParser('./data', node_name)


def rename_tree_leaves(tree):
    for cl in tree.get_terminals():
        cl.name = afn.acc2name(acc2name_path)[cl.name]
    return tree


def get_linked_HOG(prev_selected_HOG, dop_obj):
    """Get list of HOGs in newly selected node that are in common orthogroup to previous HOG"""
    in_level = prev_selected_HOG.split('.')[0]
    in_HOG_path = os.path.join(f'./data/Phylogenetic_Hierarchical_Orthogroups/{in_level}.tsv')
    in_df = pd.read_csv(in_HOG_path, sep='\t', usecols = ['HOG', 'OG']).set_index('HOG')
    orthogroup = in_df.loc[prev_selected_HOG, 'OG']
    new_HOG_list = dop_obj.OG_HOG_dict()[orthogroup]
    return new_HOG_list


def copynumber_bargraph_data_dict(HOGs, dop_obj, id_type='accession'):
    """Make a nested dict 'data_dict' {HOG{leaf_id:counts}} for input into create_tree_w_bargraphs function
    
    Parameters:
    HOGs (list): List of one or more HOGs
    dop_obj (dash_orthoparser oject): must be generated atthe same Hierarchical level as the input HOG
    id_type (str): Can be 'accession' or 'name'. Determines whether the dict keys will be accessions or sci. names.
    """
    
    data_dict =  {}
    data_dict_names = {}
    for HOG in HOGs:
        data_dict[HOG] = dop_obj.HOG_proteins_in_genome(HOG, dop_obj.accessions).to_dict()
    if id_type =='accession':
        return data_dict

    if id_type == 'name':
        for HOG in HOGs:
            data_dict_names[HOG] = {}
            for  acc, count in data_dict[HOG].items():
                name = dop_obj.accession_to_name[acc]
                data_dict_names[HOG][name] = data_dict[HOG][acc]
        return data_dict_names


def check_HOGs(HOGs, dop_obj):
    """Check if current HOG(s) are at proper 'N#' level, if not return proper HOGs from common orthogroup """
    for HOG in HOGs:
        if HOG not in dop_obj.HOGs: 
            HOGs = get_linked_HOG(HOG, dop_obj)
            if len(HOGs) > 8:
                HOGs = HOGs[:8]
            break
    return HOGs



# def get_hover_text(tree, rank_to_return=None):
#     """return a list of f strings for plotly hoverdata that correspond to each node on the tree.
    
#     Note: Terminal nodes will not be annotated beyond what is input in the clade.name field.
    
#     rank_to_return (str): Available ranks: ['no rank', 'superkingdom','family', 'genus', 'phylum', 'class', 'order', 'species',]
#     """
#     available_ranks = ['no rank', 'superkingdom','family', 'genus', 'phylum', 'class', 'order', 'species',None]
#     if rank_to_return not in available_ranks:
#         raise ValueError(f'"{rank_to_return}" not in available ranks: {available_ranks}')
        
#     with open('./data/summary_data/internal_node_dict.json', 'r') as f:
#         int_node_dict = json.load(f)
#     my_tree_clades = jdt.get_x_coordinates(tree).keys()
#     text=[]
#     for cl in my_tree_clades:
#         if cl.is_terminal():
#             text.append(cl.name)
#         else:
#             text_dict = int_node_dict.get(cl.name)
#             taxid = text_dict["taxid"]
#             sci_name = text_dict["sci_name"]
#             if rank_to_return:
#                 rank = rank_to_return
#                 rank_taxid = jnt.get_lineage_rank_dict(taxid).get(rank_to_return)
#                 if rank_taxid:
#                     sci_name = list(ncbi_tax.get_taxid_translator([rank_taxid]).values())[0]
#                 else: sci_name = None
#             else:
#                 rank = text_dict["rank"]
#             data =  (f'name: {cl.name}<br>'
#                     f'sci_name: {sci_name}<br>'
#                     f'rank: {rank}<br>'  

#                     f'taxid: {taxid}<br>'
#                     )
#             text.append(data)
#     return text

def get_hover_text(tree):
    """return a list of f strings for plotly hoverdata that correspond to each node on the tree.
    
    Note: Terminal nodes will not be annotated beyond what is input in the clade.name field.

    """
        
    with open('./data/summary_data/internal_node_dict.json', 'r') as f:
        int_node_dict = json.load(f)
    text=[]
    my_clades = tree.depths().keys()
    for node in my_clades:
        if node.is_terminal():
            text.append(node.name)
        else:
            text_dict = int_node_dict.get(node.name)
            taxid = text_dict["taxid"]
            sci_name = text_dict["sci_name"]
            rank = text_dict["rank"]
            data =  (f'name: {node.name}<br>'
                    f'sci_name: {sci_name}<br>'
                    f'rank: {rank}<br>'        
                    f'taxid: {taxid}<br>'
                )
            text.append(data)
    return text


################################################Layout###########################################################################################



app.layout = html.Div([
    ae.make_header('./assets/logo_simple.png', lab_name),
    html.H3('Empty header row'),
    ae.spTree_selector_bar(acc2name_path),
    html.Div([ # container holding keyword search and HOGs dropdown
        ae.keyword_search(acc2name_path),
        html.H5('Select HOGs to display'),
        dcc.Dropdown(id = 'HOGs_to_display_drpdwn',
                    multi=True,
        ),
        dcc.Dropdown(id = 'available_ranksdrpdwn',
                     options = [{'label':rank,'value':rank} for rank in ['superkingdom','phylum','class','order','family','genus', ]]
        ),
        dcc.Dropdown(id = 'available_rank_scinames'),


    ], style = {'width': '16%',
                'display':'inline-block',
                'align':'top',
                'float':'left',
                'padding-left':50}
    ),
    html.Div([ # container holding species tree graph
        dcc.Graph(id = 'copynumber_bargraph'),
        dcc.Store(id='node_store'),
    ],  style={
        'display':'inline-block',
        'width':'40%'}
    ),

    html.Div([# container holding species dropdowns for orthology finding
        html.H5('Find orthologs between two species'),
        dcc.Dropdown(id = 'species_ortho_1',
            options = [{'label':afn.acc2name(acc2name_path)[key],'value':key} for key in afn.acc2name(acc2name_path).keys()],
            value = None
        ),
        dcc.Dropdown(id = 'species_ortho_2',
            options = [{'label':afn.acc2name(acc2name_path)[key],'value':key} for key in afn.acc2name(acc2name_path).keys()],
            value = None
        ),
        html.Button('Download orthologs CSV', id='download-orthologs-btn'),
        dcc.Download(id="download-orthologs"),
    ],  style={
        'display':'inline-block',
        'width':'20%',
        'float':'right',
        'padding-right':50}
    ),


    #middle third
    html.Hr(),
    html.Div([#container holding HOG and species dropdowns, gene tree and genome map
        html.Div([#orthogroup explorer
            html.Div([#HOG selector dropdown
                html.H5('HOG'),
                dcc.Dropdown(id='HOG_drpdwn_genetree')                        
                ], style = {
                    'width': '100%',
                    'display':'block'}
                ),
            html.Div([ #dropdown organism selector
                html.H5('Species name'),
                dcc.Dropdown(id = 'name_drpdwn_genetree',
                    
                    options = [{'label':afn.acc2name(acc2name_path)[key],'value':key} for key in afn.acc2name(acc2name_path).keys()],
                    value = list(afn.acc2name(acc2name_path).keys())[0])                            
                ], style = {
                    'width': '100%', 
                    'display':'block'}
                ),
            html.Div(id = 'output_div',
                style = {
                    'padding-left':25,
                    'display':'block',
                    'width': '100%'}
                ),
            ], style = {
                'width': '20%',
                'padding-left':25,
                'display':'inline-block'}
            ),
        html.Div([ #gene tree and its title
            html.Div(id='title_div', 
                style={
                    'display':'block',
                    'margin':0,
                    'padding':0}
                ),     
            dcc.Graph(id = 'gene_tree_graph', 
                style={
                    'display':'block',
                    'margin':0,
                    'padding':0}
                )
        ], style={
            'display':'inline-block',
            'width': '40%',
            'vertical-align': 'top'}
        ),
        html.Div(id='stacked_genome_graphs',
            style={
                'display':'inline-block',
                'width': '35%',
                'vertical-align': 'top'}
            ),
        ], style={
            'width':'100%',
            'min-width':1500}
        ),
    #bottom third
    html.Hr(),
    html.P(),
    html.Div(id='test_output'),
])

###########################################################################################################################
 ###############################################callbacks############################################################



@app.callback(
        Output('available_rank_scinames', 'options'),
        [Input('available_ranksdrpdwn', 'value'),
        ])
def update_rank_sciname_drpdwn(rank):
    if rank:
        with open('./ranks_mrca_clades.json', 'r') as f:
            rank_mrca_clades  = json.load(f)
        mrca_clades = rank_mrca_clades[rank]
        return [{'label':sci_name,'value':cl_name} for sci_name, cl_name in mrca_clades.items()]
    else:
        return [{'label':'Select a rank','value':None}]




# @app.callback(
#         Output('output_div', 'children'),
#         [Input('name_drpdwn_genetree', 'value'),
#         Input('HOG_drpdwn_genetree', 'value'),
#         ])
# def add_uniprot_link(species_accession, HOG):
#     if HOG and species_accession:
#         tax_level=HOG.split('.')[0]
#         dop_obj = dop.DashOrthoParser('./data', tax_level=tax_level)
#         # print(dop_obj.HOG_node)
#         proteins = dop_obj.HOG_protnames_in_genome(HOG, species_accession)
#         prot_uniprot_d = pd.read_csv('proteinID_uniprotID_map.txt', 
#                     sep='\t', header=None).set_index(0).squeeze().to_dict()
#         # for prot in proteins:
#         #     print(prot, prot_uniprot_d.get(prot), species_accession)

#         protein_urls = {prot:f'https://www.uniprot.org/uniprotkb/{prot_uniprot_d.get(prot)}/entry' 
#                                                             if prot_uniprot_d.get(prot) else None for prot in proteins}
#         return html.Div([html.Li(html.A(f'uniprot link: {prot}',href=protein_urls[prot])) for prot in proteins], style={'text-align':'left'})


@app.callback(Output('download-orthologs', 'data'),
              [Input('species_ortho_1', 'value'),
               Input('species_ortho_2', 'value'),
               Input('download-orthologs-btn', 'n_clicks')])
def download_orthologs(species_1, species_2, n_clicks):
    """Download csv of orthologs between two species on button click"""
    if n_clicks:
        if species_1 != species_2:
            df = get_orthologs(species_1, species_2)
            df = get_orthology_counts(df)
            return dcc.send_data_frame(df.to_csv, filename=f'{species_1}__v__{species_2}.csv')
        else:
            dop_obj= dop.DashOrthoParser('./data')
            df = dop_obj.make_genome_annotation_df(species_1)
            return dcc.send_data_frame(df.to_csv, filename=f'{species_1}.csv')



@app.callback(
        Output('name_drpdwn_keyword', 'label'),
        Output('name_drpdwn_keyword', 'value'),
        [Input('copynumber_bargraph', 'clickData')
        ])
def names_for_kywrd_drpdwn(tree_clickData):
    acc2name_d = afn.acc2name(acc2name_path)
    tree = jdt.read_treefile(path_to_species_tree)

    def default_return():
        accs = list(acc2name_d.keys())
        value = accs[0]
        label = acc2name_d[accs[0]]
        return label, value
    
    if not tree_clickData:
        return default_return()
    elif tree_clickData['points'][0]['text'] not in acc2name_d.values():
        return default_return()
    else:
        species_name = tree_clickData['points'][0]['text']
        name2acc = {name:acc for acc, name in acc2name_d.items()}
        label = species_name
        value = name2acc[species_name]
        return label, value



@app.callback(
        Output('HOGs_to_display_drpdwn', 'options'),
        Output('HOG_drpdwn_genetree', 'options'),
        [Input('node_store', 'data'),
        ])
def HOGs_for_drpdwn_display(node_store_data):
    node_name = node_store_data['selected_node']
    tree = jdt.read_treefile(path_to_species_tree)
    int_nodes = [node.name for node in tree.get_nonterminals()]
    if node_name in int_nodes:
        dop_obj= dop.DashOrthoParser('./data', tax_level=node_name)
    else:
        dop_obj= dop.DashOrthoParser('./data', tax_level='N0')
    options = [{'label':HOG,'value':HOG} for HOG in dop_obj.HOGs]
    return options, options



    

@app.callback(
        Output('node_store', 'data'),
        [Input('available_rank_scinames', 'value'),
        Input('copynumber_bargraph', 'clickData'),
])
def update_store_with_selected_node(rank_scinames, tree_clickData):
    """Returns either node name (e.g. N38) or accesssion if node is terminal"""
    ctx = dash.callback_context
    triggered_input = ctx.triggered[0]['prop_id']
    tree = jdt.read_treefile(path_to_species_tree)
    if triggered_input == 'copynumber_bargraph.clickData':
        leaf_names = [afn.acc2name(acc2name_path)[cl.name] for cl in tree.get_terminals()] #leaf_names as sci_names
        node_name = tree_clickData['points'][0]['text']
        if node_name not in leaf_names:
            node_name = node_name.split('<br>')[0].split(' ')[1] #This will be a sci_name, since it is coming from tree data
        else:
            node_name = afn.name2acc(acc2name_path)[node_name]
        print(f'node_name in store: {node_name}')        
        return {'selected_node':node_name}
    elif triggered_input == 'available_rank_scinames.value':
        return {'selected_node':rank_scinames}
    else:
        return {'selected_node':'N0'}
    

@app.callback(
        Output('copynumber_bargraph', 'figure'),
        Output('HOGs_to_display_drpdwn', 'value'),
        [Input('HOGs_to_display_drpdwn', 'value'),
        Input('leaf_text_display', 'value'),
        Input('name_drpdwn_keyword', 'value'),
        Input('node_store', 'data'),
        ])
def species_tree_w_bargraph(HOGs,node_label_mode, species_accession, node_store):
    node_name = node_store['selected_node'] 
    if not node_name:
        node_name = 'N0'
    species_name = afn.acc2name(acc2name_path=acc2name_path)[species_accession]
    tree = jdt.read_treefile(path_to_species_tree)
    if node_name in [cl.name for cl in tree.get_nonterminals()]:
        dop_obj = dop.DashOrthoParser('./data',tax_level=node_name)
    elif node_name == 'not monophyletic':
        node_name = 'N0'
        dop_obj = dop.DashOrthoParser('./data',tax_level='N0')
    else:
        species_name = afn.acc2name(acc2name_path=acc2name_path)[node_name]
        dop_obj = dop.DashOrthoParser('./data',tax_level='N0')
    
    node_color_dict ={node_name:'rgb(250,0,0)', species_name:'rgb(250,0,0)'}
    node_size_dict = {node_name:15, species_name:15}

    print(f'node name: {node_name}\nspecies name {species_name}\ndop_node: {dop_obj.HOG_node}\nHOGs: {HOGs}')

    if HOGs:
        leaves  = [cl.name for cl in tree.get_terminals()]
        HOGs = check_HOGs(HOGs, dop_obj)
        if node_name in leaves:
            data_dict = {}
            for HOG in HOGs:
                data_dict[HOG] = {afn.acc2name(acc2name_path)[acc]:0 for acc in leaves} 
        else:
            data_dict = copynumber_bargraph_data_dict(HOGs, dop_obj,id_type='name') #name
        
        tree = rename_tree_leaves(tree) #changes tree leaves from accession to sci_name 
        text = get_hover_text(tree)
        cl_to_highlight = tree.find_any(node_name) 
        #print(f'clade_to_highlight: {node_name}')
        print(HOGs)
        return jdt.create_tree_w_bargraphs(tree, data_dict, colors, intern_node_label='name', node_color_dict=node_color_dict,cl_to_highlight=cl_to_highlight,
                                            highlight_line_width=3, label_mode=node_label_mode, hover_text=text), HOGs
    else:
        tree = rename_tree_leaves(tree)
        text = get_hover_text(tree)
        cl_to_highlight = tree.find_any(node_name)
        return jdt.create_tree(tree,  intern_node_label='name',node_color_dict=node_color_dict,node_size_dict=node_size_dict,
                               label_mode=node_label_mode, cl_to_highlight=cl_to_highlight, hover_text=text), HOGs








@app.callback(
        Output('gene_tree_graph', 'figure'),
        [Input('HOG_drpdwn_genetree', 'value'),
         Input('name_drpdwn_genetree', 'value'),
         Input('copynumber_bargraph', 'clickData'),
        ])     
def update_HOGgene_tree_graph(selected_HOG, selected_name, tree_clickData):
    ingroup = []
    outgroup = []  ## need to get rid of this and line above and args to create_HOGgene_tree()
    dop_obj = make_dopObj_from_click(tree_clickData)
    if selected_HOG:
        if selected_HOG in dop_obj.HOGs:
            orthogroup = dop_obj.HOG_OG_dict()[selected_HOG]
            path_to_gene_tree = os.path.join(path_to_gene_trees, orthogroup+'_tree.txt')
            pruned_tree = jdt.prune_tree(path_to_gene_tree, dop_obj.all_prots_in_HOG(selected_HOG))
            
            return afn.create_HOGgene_tree(pruned_tree, selected_name, ingroup, outgroup, 
                                            orthoparser_obj=dop_obj, colors=colors)
            
        else:
            trace = go.Scatter(x = list(range(100)), y = list(range(100)))
            return go.Figure(data=[trace])
    else:
        trace = go.Scatter(x = list(range(100)), y = list(range(100)))
        return go.Figure(data=[trace])     




@app.callback(
     Output('stacked_genome_graphs', 'children'),
    [Input('HOG_drpdwn_genetree','value'),
     Input('name_drpdwn_genetree','value'),
      Input('gene_tree_graph', 'clickData'),
        Input('copynumber_bargraph', 'clickData'),])    
def update_genome_map(HOG, species, tree_clickData, node_clickData ): 
    'return data, layout for plotly graph, for genes in a local area of genome '
    dop_obj = make_dopObj_from_click(node_clickData)
    if HOG:
        if HOG in dop_obj.HOGs:
            lst_of_dcc_graphs = []
            feature_type = 'CDS'
            click_protein = None
            if tree_clickData:
                click_protein = afn.get_protein_from_genetree_click(tree_clickData)
            if HOG and species:
                list_of_proteins = dop_obj.HOG_gene_dict(species)[HOG]
                for protein in list_of_proteins:
                    protein  = protein.strip()
                    fixed_prot_name = feature_type.lower()+'-'+protein
                    trimmed_df = afn.get_df_for_feature(dop_obj, species, fixed_prot_name, fts=15)
                    if protein == click_protein:
                        fig = afn.make_map(trimmed_df, fixed_prot_name, colors=colors, bgcolor='rgba(0,102,153,0.25)')
                        dcc_graph = dcc.Graph(figure=fig)
                    else:
                        dcc_graph = dcc.Graph(figure=afn.make_map(trimmed_df, fixed_prot_name, colors=colors))
                    lst_of_dcc_graphs.append(dcc_graph)   
                return html.Div(lst_of_dcc_graphs, style={'display':'block'})
        else:
            return ''
    else:
        return ''




@app.callback(
    Output('search_results', 'children'),
    [Input('name_drpdwn_keyword', 'value'),
     Input('begin_search_button','n_clicks'),
     Input('copynumber_bargraph', 'clickData'),
     State('keyword_search', 'value')
     ])
def find_protein_by_keyword(accession, n_clicks, tree_click, keyword):
    import re
    dop_obj = make_dopObj_from_click(tree_click)
    
    if n_clicks:
        df = dop_obj.make_genome_annotation_df(accession,get_common_names=True)
        pattern  = re.compile(keyword, re.IGNORECASE)
        matches = []
        for ind in df.index:
            str_to_search = str_to_search = f'{ind[4:]}, {df.loc[ind,"Parents"]}, {df.loc[ind,"Common_names"]},       \
                  {df.loc[ind,"strand"]}, {df.loc[ind,"products"]}, {df.loc[ind,"Orthogroups"]}, {df.loc[ind,"HOGs"]}'
            if pattern.search(str_to_search):
                matches.append(str_to_search)
        if not matches:
            return html.H5(f'There were no matches for {keyword} in {dop_obj.accession_to_name[accession]}')
        else:
            return html.Div([html.Li(ele) for ele in matches], style={'text-align':'left'})
    else:
        return html.H4(f'pick organism and search its protein product description with a keyword...')


if __name__ == '__main__':
    app.run_server()
