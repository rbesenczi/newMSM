import pandas as pd
import plotly.express as px
from os import environ

home = environ['HOME']
dataset = "HCP"
data = pd.read_csv(home + '/groupwise/' + dataset + '/stats.csv')

fig = px.scatter(data, x='subject', y=['typical','group'], title="Individual distortion")
fig.update_layout(
    xaxis_title="Subjects",
    yaxis_title="Areal 95% distortion",
    legend_title="Type of registration"
)

fig.update_traces(textposition="bottom right")

fig.show()
