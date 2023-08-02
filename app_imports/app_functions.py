from Bio import Phylo
import pandas as pd
import os
import json
import plotly.graph_objs as go
import numpy as np


def calc_wink_value(your_HOG_enrich, all_HOG_enrighments):
    import bisect
    all_HOG_enrighments_lst =list(all_HOG_enrighments.sort_values('Enrichment')['Enrichment'])
    index = bisect.bisect(all_HOG_enrighments_lst, your_HOG_enrich)
    p_t = index/all_HOG_enrighments.shape[0]
    if p_t>=0.5:
        p =1-p_t
    else:
        p=p_t
    p = (int(p*1000))/1000
    return  p



def Enrichment_histogram(ingroup, outgroup, copies, HOG, dop_obj):
    enrich_dict = dop_obj.HOG_enrichment(ingroup,outgroup, copies=copies)
    df = pd.DataFrame().from_dict(enrich_dict, orient='index')
    df.columns = ['Enrichment']
    df = df.sort_values('Enrichment')
    counts, bins = np.histogram(df['Enrichment'], bins=range(-105,105, 5), density=False)
    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=bins, 
        y=counts,
        name = 'Enrich. distribution'))
    
    if HOG in enrich_dict.keys():
        val = enrich_dict[HOG]
        wink_value = calc_wink_value(val, df)
        fig.add_trace(go.Scatter(
            x = [val],
            y = [0],
            name = f'{HOG}<br>p={wink_value}'
            ))
    fig.update_layout(
        paper_bgcolor='rgb(250,250,250)',   
        plot_bgcolor='rgb(253,253,253)',
        margin={'t':20, 'b':0, 'l':0, 'r':0}
        )
    fig.update_xaxes(
        title_text = "HOG enrichment",
        title_font = {"size": 12},
        tickfont = {'size':8}
        )
    fig.update_yaxes(
        title_text = "Counts",
        title_font = {"size": 12},
        tickfont = {'size':8}
        )
    
    return fig



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
    


def filter_enrichment_values(enrichment_dict, threshold):
    #keep enrichment # above or below threshold
    filtered_enrichment_dict = {}
    min_cut = -threshold
    max_cut = threshold
    for key in enrichment_dict.keys():
        is_between = min_cut <= enrichment_dict[key] <= max_cut
        if not is_between:
            filtered_enrichment_dict[key] = enrichment_dict[key]
    return filtered_enrichment_dict



def h_or_o_enrichments(ingroup, outgroup, copies, threshold, annot_species = None, group_type = 'HOG', orthoparser_obj=None):
    if annot_species:
        species_to_pull_annots = annot_species
    elif len(ingroup)>0:
        species_to_pull_annots = ingroup[0]
    else:
        annot_species = None
        
    if group_type == 'HOG':
        enrichments = orthoparser_obj.HOG_enrichment(ingroup,outgroup,copies)
        
    filtered_enrichment = filter_enrichment_values(enrichments, threshold)
    df = pd.DataFrame(list(filtered_enrichment.values()), list(filtered_enrichment.keys()) )
    df.columns = ['Enrichment']
    df = df.sort_values('Enrichment',ascending = False)   
    df_t = pd.read_csv(f'./data/genome_annotations/{species_to_pull_annots}_annotations.csv')
    df_t = df_t.set_index('HOGs')
    HOGs = list(df.index)
    functions = []
    for HOG in HOGs:
        if HOG in df_t.index:
            function = df_t.loc[HOG,'products']
            if not type(function) ==str:
                functions.append(function[0])  
            else:
                functions.append(function)
        else:
            functions.append('nan')
    return df, functions



def create_HOGgene_tree(tree, selected_name, ingroup, outgroup, orthoparser_obj=None, colors=None):

    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
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



def HOG_enrichment_figure(ingroup,outgroup,copies, threshold, annot_species, orthoparser_obj=None, colors=None):
    #df, functions =  HOG_enrichment(ingroup, outgroup, copies, threshold, annot_species)
    df, functions = h_or_o_enrichments(ingroup, outgroup, copies, threshold, annot_species, group_type = 'HOG', orthoparser_obj=orthoparser_obj)
    trace1 =go.Bar(x = df.index, 
                   y = df['Enrichment'],
                   text = functions,
                   hoverinfo = 'text',
                   marker_color = colors['t_green'],
                   #marker_line_color = 'rgb(240,240,240)',
                   #marker_line_width=1.5,
                   opacity=0.6
                   )
    data = [trace1]
    layout = go.Layout({'paper_bgcolor':'rgb(250,250,250)',
                        'plot_bgcolor':'rgb(250,250,250)',
                        'margin':{'b': 0, 'l': 60, 'r': 10, 't': 0}})          
    fig = go.Figure(data, layout)
    fig.update_traces(marker_color=colors['t_blue'],) #marker_line_color='grey',
                  #marker_line_width=2, opacity=0.8)
    
    return fig      



def create_tree(tree_filepath, ingroup, outgroup, HOG, copies, title=None, orthoparser_obj=None, colors=None):
    missing_dict = orthoparser_obj.genomes_missing_HOGs(ingroup, outgroup, HOG, copies)
    tree = read_treefile(tree_filepath)
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    name_dict = orthoparser_obj.accession_to_name
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []
    node_sizes = []
    color_dict = {}
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        if cl.name in name_dict.keys():
            text.append(name_dict[cl.name])
        else:
            text.append(cl.name)
            
        if cl.name in ingroup:
            color_dict[cl.name] = colors['t_green']  #seagreen trestle
            node_sizes.append(10)
        elif cl.name in outgroup:
            color_dict[cl.name] = colors['t_blue']
            node_sizes.append(10)
        else:
            color_dict[cl.name] = '#92A0A9' #grey'
            node_sizes.append(5)
        if cl.name in missing_dict['ingroup']:
           color_dict[cl.name] = 'rgba(150,255,150,0.65)'
        if cl.name in missing_dict['outgroup']:
           color_dict[cl.name] = 'rgba(0,130,255,0.4)'

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
                marker=dict(color=list(color_dict.values()),
                            size=node_sizes),
                text=text,  # vignet information of each node
                hoverinfo='text',
                )
    if title:
        title=title
    layout = dict(title=title,
                  paper_bgcolor='rgb(248,248,248)',
                  dragmode="lasso",
                  font=dict(family='Balto', size=14),
                  #width=750,
                  height=550,
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
                  plot_bgcolor='rgb(248,248,248)',
                  legend={'x': 0, 'y': 1},
                  margin={'b': 0, 'l': 0, 'r': 0, 't': 0}
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



def get_df_for_feature(assembly_accession, feature_type, feature):
    'Return a local df containing features near the passed feature'
    
    filename = assembly_accession +  '_annotations.csv'
    df = pd.read_csv(f'./data/genome_annotations/{filename}') 
    df = df.set_index('IDs')
    fts=15 # fts = features to show on each side of of feature of interest
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
    


def prune_tree(path_to_tree, leaves_to_keep):
    tree = read_treefile(path_to_tree)
    pruned_tree = read_treefile(path_to_tree)
    for leaf in tree.get_terminals():
        if leaf.name not in leaves_to_keep:
            pruned_tree.prune(leaf.name)
    leaves_removed= tree.count_terminals() - pruned_tree.count_terminals()
    return pruned_tree



def get_accession_from_genetree_click(geneTree_clickData, dop_obj):
    name = ' '.join(geneTree_clickData['points'][0]['text'].split(' ')[:-1]).strip() 
    return dop_obj.name_to_accession[name]


def get_protein_from_genetree_click(geneTree_clickData):
    return geneTree_clickData['points'][0]['text'].split(' ')[-1]
################################################################################################################################################################
#
#TREE DRAWING
#
################################################################################################################################################################
def get_x_coordinates(tree):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth
    # is the distance from root to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords


def get_y_coordinates(tree, dist=1.3):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords


def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=0.5):
    """define a shape of type 'line', for branch
    """
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")

    return branch_line


def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=1, x_coords=0, y_coords=0):
    """Recursively draw the tree branches, down from the given clade"""

    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)

    line_shapes.append(branch_line)

    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]

        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)


def read_treefile(filename):
    'create tree object from newick format using Bio.Phylo'
    tree = Phylo.read(filename, "newick")
    return tree

################################################################################################################################################################
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#TREE DRAWING
#
################################################################################################################################################################
def create_simple_tree(tree_filepath, title='Tree'):
    tree = read_treefile(tree_filepath)
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords,
               y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []
    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
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
                #marker=dict(color=list(color_dict.values()),
                            #size=node_sizes),
                #text=text,  # vignet information of each node
                #hoverinfo='text',
                )
    layout = dict(title=title,
                  paper_bgcolor='rgb(248,248,248)',
                  dragmode="lasso",
                  font=dict(family='Balto', size=14),
                  #width=750,
                  height=550,
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
                  plot_bgcolor='rgb(248,248,248)',
                  legend={'x': 0, 'y': 1},
                  margin={'b': 0, 'l': 0, 'r': 0, 't': 0}
                  )
    fig = dict(data=[data], layout=layout)
    return fig 