
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

pio.renderers.default = "browser"

# Plane
x = np.linspace(0, 2, 1000)
y = np.linspace(0, 2, 1000)
X, Y = np.meshgrid(x, y)
Z = 1 - X - Y
Z[X + Y > 1] = np.nan

fig = go.Figure()

# Plane
fig.add_trace(go.Surface(x=X, y=Y, z=Z, opacity=0.5, showscale=False, colorscale=[[0, '#00fa47'], [1, '#00fa47']],
                         contours={ "x": {"show": True, "start":X.min(), "end":X.max(), "color": "#00fa47", "size": 0.05},
                                    "y": {"show": True, "start":Y.min(), "end":Y.max(), "color": "#00fa47", "size": 0.05},
                                    "z": {"show": True, "start":Z.min(), "end":Z.max(), "color": "#00fa47", "size": 0.02}}))

x1 = np.linspace(0, 2, 1000)
z1 = np.linspace(0, 2, 1000)
X1, Z1 = np.meshgrid(x1, z1)

k1 = 0.8
k2 = 1.2
a1 = k1/(k1+k2)
a2 = k2/(k1+k2)
L = 0.1
Y1 = (a1*X + L)/a2
Z1[(X1>0.6) | (Y1>0.6)] = np.nan

fig.add_trace(go.Surface(x=X1, y=Y1, z=Z1, opacity=0.5, showscale=False, colorscale=[[0, '#d400c6'], [1, '#d400c6']],
                         contours={ "x": {"show": True, "start":X1.min(), "end":X1.max(), "color": "#d400c6", "size": 0.05},
                                    "y": {"show": True, "start":Y1.min(), "end":Y1.max(), "color": "#d400c6", "size": 0.05},
                                    "z": {"show": True, "start":Z1.min(), "end":Z1.max(), "color": "#d400c6", "size": 0.02}}))

# x-axis
fig.add_trace(go.Scatter3d(x=[0, 1.3], y=[0, 0], z=[0, 0], mode="lines", line=dict(color="white", width=6), name="RONO2"))

# y-axis
fig.add_trace(go.Scatter3d(x=[0, 0], y=[0, 1.3], z=[0, 0], mode="lines", line=dict(color="white", width=6), name="NO2"))

# z-axis
fig.add_trace(go.Scatter3d(x=[0, 0], y=[0, 0], z=[0, 1.3], mode="lines", line=dict(color="white", width=6), name="NO"))


t = np.linspace(0, 5/6, 100000) 
X2 = 0.5 - 3*t/5
Y2 = 0.5 - 2*t/5
Z2 = t

fig.add_trace(go.Scatter3d(
    x=X2, y=Y2, z=Z2,
    mode="lines",
    line=dict(width=8, color='white'),
    name="intersection line",
    hoverinfo='skip'
))

fig.add_trace(go.Scatter3d(
    x=[1.2, 0, 0],
    y=[0, 1.2, 0],
    z=[0, 0, 1.2],
    mode="text",
    text=["<b>RONO2</b>", "<b>NO2</b>", "<b>NO</b>"],
    textfont=dict(color="white", size=14),
    showlegend=False
))


fig.update_layout(scene=dict(bgcolor='black',
                xaxis=dict(title="", range=[0,1.2], showgrid=False, zeroline=False, showbackground=False, showticklabels=False),
                yaxis=dict(title="", range=[0,1.2], showgrid=False, zeroline=False, showbackground=False, showticklabels=False),
                zaxis=dict(title="", range=[0,1.2], showgrid=False, zeroline=False, showbackground=False, showticklabels=False), 
                annotations=[
                dict(
                        showarrow=True,
                        x=0.28,   
                        y=0.34,  
                        z=0.4,    
                        text="Intersection Line",
                        textangle=0,
                        ax=100,  
                        ay=-100,
                        arrowhead=2,
                        arrowcolor="red",
                        arrowsize=1,
                        arrowwidth=2,
                        font=dict(color="black", size=12),
                        bgcolor="yellow",
                        bordercolor="black",
                        borderwidth=1)]),
                )


fig.show()