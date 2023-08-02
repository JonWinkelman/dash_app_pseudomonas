import dash
from dash import dcc, html
import plotly.graph_objs as go
import base64
from app_imports import app_functions_dev as afn


def make_header(image_filepath, lab_name):

    encoded_image = base64.b64encode(open(image_filepath, 'rb').read())
    image = html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()),
                height=100, width=100, id='logo',
                style={
                    'padding-top':0,
                    'padding-left':50,
                    'display':'inline-block'}
                )

    return html.Div([
                html.Div([
                    image
                    ],  
                    className='col-1',
                    style={
                        'width': '33%',
                        'display':'inline-block',
                        'float':'left'}
                    ),
                html.Div([
                    html.H1(['Hierarchical Orthogroups'], id='page-title')
                    ],  
                    className="col-2",
                    style={
                        'width': '33%',
                        'display':'inline-block',
                        'vertical-align':'top'}
                    ),
                html.Div([
                    html.H1([lab_name], 
                        id='lab_name',
                        style={ 
                            'text-align':'right',
                            'padding-right':100}
                        ),
                    ], 
                    className="col-3",
                    style={
                        'width': '33%',
                        'display':'inline-block',
                        'float':'right',
                        'vertical-align':'top',
                        'min-width':300}
                    ),
                ], 
                className="header-container",
                style={
                    'width': '100%',
                    'padding-top':15,
                    'padding-bottom':0,
                    'overflow':'hidden',
                    'background-color':'rgb(248,248,248)',
                    'min-width':900}
                )



                         
def spTree_selector_bar(acc2name_path):
    leaf_label_button = html.Div([
    # html.Div([
    #     dcc.RadioItems(id = 'reset_button',
    #         options=[
    #             {'label': 'original', 'value': 'original'},
    #             {'label': 'prune', 'value': 'prune'}],
    #             value='original',
    #             labelStyle={'display':'flex'}),
    #     ], style = {          
    #         'display': 'inline-block',
    #         'width' :'30%',
    #         'padding':20}
    #     ),
    html.Div([
        dcc.RadioItems(id = 'leaf_text_display',
            options=[
                {'label': 'No labels', 'value': 'markers'},
                {'label': 'Labels', 'value': 'markers+text'}],
                value='markers',
                labelStyle={'display':'flex'}),

        ], style = {
            'display': 'inline-block',
            'width' :'30%',
            'padding':20}
        ),
    # html.Div([

    #     dcc.Dropdown(id = 'species_to_keep',
    #         options = [{'label':afn.acc2name(acc2name_path)[key],'value':key} for key in afn.acc2name(acc2name_path).keys()],
    #         multi=True,
    #         value = list(afn.acc2name(acc2name_path).keys())[:2]),
    #     ], style = {
    #         'display': 'inline-block',
    #         'width' :'30%',
    #         'padding':20}
    #     ),

    ])
    return leaf_label_button

              



def keyword_search(acc2name_path):

    keyword_search = html.Div([
                html.H4('keyword search'),
                dcc.Dropdown(id = 'name_drpdwn_keyword',
                     options = [{'label':afn.acc2name(acc2name_path)[key],'value':key} for key in afn.acc2name(acc2name_path).keys()],
                     value = list(afn.acc2name(acc2name_path).values())[0]),
                dcc.Input(
                    id='keyword_search',
                    type = 'text',
                    placeholder='type search term here...',
                    style = {'width': '100%', 'display': 'inline-block'}
                    ),
                html.Button("GO", id="begin_search_button"),
                html.Div(id='search_results')
                    ],style={ 
                        'display': 'inline-block',
                        'width':'100%',
                        }
                    )
    return keyword_search
                
