import dash
from dash import dcc, html
import plotly.graph_objs as go
import base64


def make_header(image_filepath, lab_name):

    encoded_image = base64.b64encode(open(image_filepath, 'rb').read())
    image = html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                    height=100, width=100, id='logo',
                    style={'padding-top':0,
                            'padding-left':50,
                            'display':'inline-block',
                            })

    return html.Div([
                html.Div([
                    image
                ],  className='col-1',
                    style={
                        'width': '33%',
                        'display':'inline-block',
                        'float':'left',
                        }),
                html.Div([
                    html.H1(['Orthogroups'], id='page-title')
                ],  className="col-2",
                    style={
                        'width': '33%',
                        'display':'inline-block',
                        'vertical-align':'top',
                    }
                ),
                html.Div([
                html.H1([lab_name], 
                    id='lab_name',
                    style={ 
                        'text-align':'right',
                        'padding-right':100,      
                    }
                ),
                ], className="col-3",
                    style={
                        'width': '33%',
                        'display':'inline-block',
                        'float':'right',
                        'vertical-align':'top',
                        'min-width':300
                        }),
        ], className="header-container",
            style={
                'width': '100%',
                'padding-top':15,
                'padding-bottom':0,
                'overflow':'hidden',
                'background-color':'rgb(248,248,248)',
                'min-width':900})




def ingroup(dop_obj, colors):
    ingroup = html.Div([
                    html.H4('Ingroup species'),
                    dcc.Dropdown(id = 'ingroup_drpdwn',
                                options = [{'label':dop_obj.accession_to_name[key],'value':key} for key in dop_obj.accession_to_name.keys()],
                                multi=True,
                                value = list(dop_obj.accession_to_name.keys())[:2]),
                    ], style = {'width': '46%',
                                'display':'inline-block',
                                'padding':15,
                                'border-style': 'solid',
                                'border-radius': 18,
                                'border-color': colors['t_green'],
                                'margin':15,
                                'margin-right':7.5,
                                'vertical-align':'top',
                                })
    return ingroup



def outgroup(dop_obj, colors):
    """make dropdown with outgroup species to select"""
    outgroup = html.Div([
                    html.H4('Outgroup species'), 
                    dcc.Dropdown(id = 'outgroup_drpdwn',
                                options = [{'label':dop_obj.accession_to_name[key],'value':key} for key in dop_obj.accession_to_name.keys()],
                                multi=True, 
                                value = list(dop_obj.accession_to_name.keys())[-2:],
                    )], style = {'width': '46%',
                                'display':'inline-block',
                                'padding':15,
                                'border-style': 'solid',
                                'border-radius': 18,
                                'border-color': colors['t_blue'],
                                'margin':15,
                                'margin-left':7.5,
                                'display':'inline-block',
                                #'float':'right',
                                'vertical-align':'top',})
    return outgroup




def sliders_drpdwn(dop_obj):
    sliders_drpdwn = html.Div([#container holding sliders, dropdown and associated titles
            html.Div([
                    html.H5('Select number of gene copies'),
                    dcc.Slider(id='copy_slider',
                                min = 0, 
                                max = 10,
                                step = 1,
                                value = 2,
                                marks = {num:num for num in range(11)}),
                    ], style = {'width': '100%',
                                'display':'block',
                                'padding-left':0,}),
            html.Div([
                html.H5('Set threshold for display'),
                dcc.Slider(id='bargraph_threshold',
                            min = 1, 
                            max = 99,
                            step = 1,
                            value = 10,
                            marks = {(num*10):(num*10) for num in range(1, 10)})
                    ], style = {'width': '100%',
                                'display':'block',
                                'padding-left':0,
                                }),
                html.Div([
                html.H5('Select genome to pull annotations from'),
                dcc.Dropdown(id='genome_for_annotations',
                                options = [{'label':dop_obj.accession_to_name[key],'value':key} for key in dop_obj.accession_to_name.keys()],
                                multi = False,
                                value = list(dop_obj.accession_to_name.keys())[0]
                            )
                    ], style = {'width': '100%',
                                'display':'block',
                                'padding-left':0}),
                html.Div([      
                            html.Button("Download data", id="btn-download-txt"), #styled in external css file
                            dcc.Download(id="download-ingrp-csv"),
                            dcc.Download(id="download-outgrp-csv"),
                            dcc.Download(id='download-species_acc-csv'),
                        ], style = {'width': '100%',
                                    'display':'block',
                                    'margin-left':'auto',
                                    'margin-right':'auto'}),
                            
            ], style = {'width': '16%',
                        'padding': 0,
                        'padding-top':120,
                        'min-width':200,
                        'display':'inline-block',
                        #'border':'solid'
                        })
    return sliders_drpdwn



def enrichment_graph():
    enrichment_graph = html.Div([
                            html.H4(['Percentage enrichment/depletion'], style={'text-align':'center'}),
                            dcc.Graph(id = 'Enrichment_graph',
                            style = {'width': '100%',
                                    'padding-top': 20}),
                        ], style = {'width': '45%',
                                    'padding': 0,
                                    'display':'inline-block',
                                    'min-width':300,
                                    #'border':'solid',
                                    'vertical-align':'top'})
    return enrichment_graph



def genome_graph():
    genome_graph = html.Div([
                        dcc.Graph(id = 'genome_graph', style={'display':'block'}),
                        dcc.Graph(id = 'genome_graph2', style={'display':'block'}),
                        dcc.Graph(id = 'genome_graph3', style={'display':'block'}),
                    ], style={
                        'padding-right':0,
                        'padding-left': 0})
    return genome_graph



def orthgrp_fndr(dop_obj):

    orthgrp_fndr = html.Div(
                        [
                            html.Div([
                            html.H4('Organism'),
                            dcc.Dropdown(
                                id='species_selector',
                                options=[{'label':dop_obj.accession_to_name[key], 'value':key} for key in dop_obj.accession_to_name.keys()],
                                value = list(dop_obj.accession_to_name.keys())[0]
                                ),
                                ],style={'width': '30%', 'display': 'inline-block'}),
                            
                            html.Div([
                            html.H4('Gene/protein list'),      
                            dcc.Dropdown(
                                id='opt_feature_dropdown',
                                ),
                                ],style={'width': '20%', 'display': 'inline-block'}
                            ),
                            html.P(),
                            html.Div(id='display_selected_values'),
                            html.Hr(),
                        ]
                    )
    return orthgrp_fndr                 


def keyword_search():

    keyword_search = html.Div([
                html.H4('keyword search'),
                dcc.Input(
                    id='keyword_search',
                    type = 'text',
                    placeholder='type search term here...',
                    style = {'width': '30%', 'display': 'inline-block'}
                    ),
                html.Button("GO", id="begin_search_button"),
                html.Div(id='search_results')
                    ],style={ 
                        'display': 'inline-block',
                        'width':'100%'
                        }
                    )
    return keyword_search
                
                         








