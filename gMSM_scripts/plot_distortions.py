import pandas as pd
import plotly.express as px
from os import environ

home = environ['HOME']
dataset = "HCP"
data = pd.read_csv(home + '/groupwise/' + dataset + '/stats.csv')

fig = px.scatter(data, x='subject', y=['global','group'])
fig.update_traces(textposition="bottom right")

fig.show()
