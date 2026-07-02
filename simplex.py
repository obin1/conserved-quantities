# import numpy as np
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure(figsize=(6,6))
# ax = fig.add_subplot(111, projection='3d')

# # Limits
# # L = 2
# # ax.set_xlim(-L, L)
# # ax.set_ylim(-L, L)
# # ax.set_zlim(-L, L)

# # Draw custom axes
# ax.quiver(0, 0, 0, 2, 0, 0, color='red', arrow_length_ratio=0.08)
# ax.quiver(0, 0, 0, 0, 2, 0, color='green', arrow_length_ratio=0.08)
# ax.quiver(0, 0, 0, 0, 0, 2, color='blue', arrow_length_ratio=0.08)

# # Label the ends
# # ax.text(L, 0, 0, 'x')
# # ax.text(0, L, 0, 'y')
# # ax.text(0, 0, L, 'z')

# x = np.linspace(0,1,100)
# y = np.linspace(0,1,100)

# X,Y = np.meshgrid(x,y)
# Z = 1 - X - Y

# ax.plot_surface(X, Y, Z, alpha=0.5)

# # plane_surface = ax.plot_surface(X, Y, Z, alpha=0.8)

# ax.xaxis.pane.fill = False
# ax.yaxis.pane.fill = False
# ax.zaxis.pane.fill = False

# ax.grid(False)
# ax.set_box_aspect([1, 1, 1])

# plt.show()

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
fig.add_trace(go.Surface(x=X, y=Y, z=Z, opacity=0.6, showscale=False, colorscale=[[0, '#00fa47'], [1, '#00fa47']]))

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

fig.add_trace(go.Surface(x=X1, y=Y1, z=Z1, opacity=0.6, showscale=False, colorscale=[[0, '#d400c6'], [1, '#d400c6']]))

# x-axis
fig.add_trace(go.Scatter3d(x=[0, 1.2], y=[0, 0], z=[0, 0], mode="lines", line=dict(color="red", width=6), name="RONO2"))

# y-axis
fig.add_trace(go.Scatter3d(x=[0, 0], y=[0, 1.2], z=[0, 0], mode="lines", line=dict(color="green", width=6), name="NO2"))

# z-axis
fig.add_trace(go.Scatter3d(x=[0, 0], y=[0, 0], z=[0, 1.2], mode="lines", line=dict(color="blue", width=6), name="NO"))


fig.update_layout(scene=dict(xaxis=dict(title="RONO2", range=[0,1.2], showgrid=False, zeroline=False, showbackground=False, showticklabels=False),
                yaxis=dict(title="NO2", range=[0,1.2], showgrid=False, zeroline=False, showbackground=False, showticklabels=False),
                zaxis=dict(title="NO", range=[0,1.2], showgrid=False, zeroline=False, showbackground=False, showticklabels=False)))


fig.show()