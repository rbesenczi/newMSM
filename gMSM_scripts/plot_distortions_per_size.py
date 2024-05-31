import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
from os import environ

home = environ['HOME']
dataset = "HCP"
data = pd.read_csv(home + '/groupwise/' + dataset + '/dist.csv')

data["Group"] = data["Group"].astype(str)
data["Lambda"] = data["Lambda"].astype(str)

fig = px.scatter(data, y='Areal 95%', x='Size', title='Distortion in function of group sizes for different lambdas', color='Lambda', symbol="Lambda", trendline="ols")
fig.update_traces(textposition="bottom right")

fig.show()
