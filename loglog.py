#%%
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
pio.renderers.default = "browser"

survey = {
    "Mechanism": ["JPMv0.2", "JPM1.1", "SmallStrato", "Superfast", "Form 1985", "POLLU", "GC-Hg", "E3SM", "CBM-Z", "SAPRC99", "CB05", "RACM", "RADM2", "MOZART", "MOZART-T1", "AMORE", "PACT-1D", "RCIM", "CIM", "CIM-var", "JAM", "GCv12.0", "GCv12.7", "GCv12.8", "GCv12.9", "GCv13.3", "GCv13.4", "GCv14.6", "CRI", "MECCA", "MCM", "Toluene", "Log81", "CRACMM2", "CRACMM3"],
    "#Species": [11, 16, 5, 15, 12, 20, 32, 47, 67, 74, 82, 82, 59, 81, 155, 133, 167,145, 386, 394, 245, 235, 243, 258, 262, 287, 287, 353, 442, 733, 5832, 12, 38, 196, 226],
    "#Reactions": [10,13, 10, 32, 25, 25, 94, 104, 142, 211, 205, 250, 156, 196, 360, 330, 513, 379, 886, 886, 702, 725, 750, 825, 850, 903, 913, 1058, 1258, 2323, 13140, 10, 57, 531, 614],
    "Stoichiometric Invariants": [2, 5, 1, 0, 1, 3, 2, 1, 7, 1, 7, 1, 2, 0, 0, 0, 9, 0, 0, 2, 4, 6, 6, 6, 8, 10, 10, 9, 2, 9, 1, 3, 10, 4, 5],
    "Coproduction Index": [1, 1, 1, 6, 6, 4, 65, 19, 10, 13, 14, 13, 9, 17, 41, 42, 49, 85, 650, 474, 53, 136, 143, 190, 188, 185, 181, 168, 216, 456, 3558, 1, 11, 22, 35],
    "Kinetic Invariants": [1, 1, 0, 1, 0, 0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 150, 17, 0, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 1, 1, 0, 0],
    "Broken Null Cycles": [0, 0, 1, 5, 6, 4, 61, 18, 10, 13, 14, 13, 9, 17, 41, 42, 48, 85, 500, 457, 53, 135, 142, 189, 186, 183, 181, 168, 216, 456, 3558, 0, 10, 22, 35]
}

mech_survey = pd.DataFrame(survey)

mech_survey["R-gamma"] = mech_survey.loc[:, "#Reactions"] - mech_survey.loc[:, "Coproduction Index"]

# remove MCM, MECCA, Toluene
mech_survey1 = mech_survey[~mech_survey["Mechanism"].isin(["MCM", "MECCA", "Toluene"])]

# Only those with KI
mech_survey2 = mech_survey1[mech_survey1["Kinetic Invariants"] != 0]

# remove CRI and GC 14.6
mech_survey3 = mech_survey1[mech_survey1["R-gamma"] < 800]


#%%
# -----------------------
# Plotly for Zooming In
# -----------------------

# Create Figure
fig = go.Figure()

# specific labels, positions, and sizes per mechanism
label_mechs = {"GCv13.4", "JAM", "CRACMM3", "CRACMM2", "MOZART-T1", "RCIM", "AMORE", "RACM", "SAPRC99", "CB05", "MOZART", "RADM2", "CBM-Z"}  
labels = [name if name in label_mechs else ""
    for name in mech_survey["Mechanism"]]
positions = ["middle right" if name in {"RCIM", "CB05", "CBM-Z"} else "middle left"
    for name in mech_survey["Mechanism"]]
sizes = [8 if name in {"Form 1985", "POLLU", "SmallStrato"}
         else 12 for name in mech_survey["Mechanism"]]

# Scatter: mech_survey3 (regular points), main plot
fig.add_trace(
    go.Scatter(
        x=np.log10(mech_survey["#Species"]),
        y=np.log10(mech_survey["R-gamma"]),
        mode="markers+text",
        marker=dict(size=sizes, color="gray", opacity=1),
        text=labels,
        textposition=positions,
        textfont=dict(size=10, color='darkgray'),
        name="Mechanisms w/o Kinetic Invariants"))

# specific labels, positions, sizes per mechanism
# label_mechs2 = {"GC-Hg", "Superfast", "JPM1.1", "JPMv0.2", "Log81"} 
# labels2 = ["" if name in label_mechs2 else name
#     for name in mech_survey2["Mechanism"]]
positions1 = ["top center" if name in {"E3SM"} 
              else "bottom center" if name in {"CIM"}
              else "middle right" for name in mech_survey2["Mechanism"]]
sizes1 = [13 if name in {"Log81", "Superfast", "GC-Hg", "JPM1.1", "JPMv0.2"}
          else 25 for name in mech_survey2["Mechanism"]]

# Scatter: mech_survey2 (stars), main plot
fig.add_trace(
    go.Scatter(
        x=np.log10(mech_survey2["#Species"]),
        y=np.log10(mech_survey2["R-gamma"]),
        mode="markers+text",
        marker=dict(
            symbol="star",
            size=sizes1,
            color="purple",
            opacity=1,
            line=dict(color="black", width=1.2)),
        text=mech_survey2["Mechanism"],
        textposition=positions1,
        textfont=dict(size=16),
        name="Mechanisms w/ Kinetic Invariants"))

# 1-1 line, main plot
x_line = np.linspace(0, np.log10(mech_survey["#Species"]).max(), 100)
y_line = x_line
fig.add_trace(
    go.Scatter(
        x=x_line,
        y=y_line,
        mode="lines",
        line=dict(color="red", dash="dash"),
        name="1–1 line",
        showlegend=False))

# Fill under 1-1 line, main plot
fig.add_trace(
    go.Scatter(
        x=np.concatenate([x_line, x_line[::-1]]),
        y=np.concatenate([np.zeros_like(x_line), y_line[::-1]]),
        fill="toself",
        fillcolor="rgba(128,0,128,0.2)",
        line=dict(color="rgba(255,255,255,0)"),
        hoverinfo="skip",
        showlegend=False))

# Layout properties, main plot
fig.update_layout(
    xaxis_title="# Species",
    yaxis_title="Effective Reactions (R-γ)",
    xaxis=dict(
        title_font=dict(size=20), 
        tickfont=dict(size=16)),
    yaxis=dict(
        title_font=dict(size=20),  
        tickfont=dict(size=16)    ),
    template="simple_white",
    width=800,
    height=800)
fig.update_xaxes(range=[0, 800])
fig.update_yaxes(range=[0, 800])

# # filter only mechanisms in crowded bottom left corner
# mech_survey4 = mech_survey[mech_survey["Mechanism"].isin(["GC-Hg", "Superfast", "JPM1.1", "JPMv0.2", "Log81"])]
# mech_survey5 = mech_survey[mech_survey["Mechanism"].isin(["SmallStrato", "Form 1985", "POLLU"])]
# fig1 = go.Figure()

# # specific positions per mechanism, inset
# positions2 = [
#     "middle right" if name in {"GC-Hg", "JPM1.1"} 
#     else "top center" if name in {"Superfast"}
#     else "middle left" if name in {"Log81"}
#     else "bottom center" 
#     for name in mech_survey4["Mechanism"]]

# # star points, inset
# fig1.add_trace(
#     go.Scatter(
#         x=mech_survey4["#Species"],
#         y=mech_survey4["R-gamma"],
#         mode="markers+text",
#         marker=dict(
#             symbol="star",
#             size=25,
#             color="purple",
#             line=dict(color="black", width=1.2)),
#         text=mech_survey4["Mechanism"],
#         textposition=positions2,
#         textfont=dict(size=16),
#         name="Mechanisms w/ Kinetic Invariants",
#         showlegend=False))

# # specific positions per mechanism, inset
# positions3 = [
#     "middle right" if name in {"POLLU"} 
#     else "middle left" if name in {"Form 1985"}
#     else "top center" 
#     for name in mech_survey5["Mechanism"]]

# # gray points, inset
# fig1.add_trace(
#     go.Scatter(
#         x=mech_survey5["#Species"],
#         y=mech_survey5["R-gamma"],
#         mode="markers+text",
#         marker=dict(size=12, color="gray"),
#         text=mech_survey5["Mechanism"],
#         textposition=positions3,
#         textfont=dict(size=8, color='darkgray'),
#         name="Mechanisms w/o Kinetic Invariants",
#         showlegend=False))

# # 1-1 line, inset
# x_line = np.linspace(0, mech_survey["#Species"].max(), 100)
# y_line = x_line
# fig1.add_trace(
#     go.Scatter(
#         x=x_line,
#         y=y_line,
#         mode="lines",
#         line=dict(color="red", dash="dash"),
#         name="1–1 line",
#         showlegend=False))

# # fill under 1-1 line, inset
# fig1.add_trace(
#     go.Scatter(
#         x=np.concatenate([x_line, x_line[::-1]]),
#         y=np.concatenate([np.zeros_like(x_line), y_line[::-1]]),
#         fill="toself",
#         fillcolor="rgba(128,0,128,0.2)",
#         line=dict(color="rgba(255,255,255,0)"),
#         hoverinfo="skip",
#         showlegend=False))

# # layout properties, inset
# fig1.update_layout(
#     xaxis=dict(tickfont=dict(size=16)),
#     yaxis=dict(tickfont=dict(size=16)),
#     template="simple_white",
#     width=800,
#     height=800,
#     showlegend=False)

# # define ranges for inset
# x0, x1 = 0, 50
# y0, y1 = 0, 50

# # Define the inset axes (xaxis2, yaxis2) inside the main figure
# fig.update_layout(
#     xaxis2=dict(
#         domain=[0.58, 0.98],   
#         anchor="y2",
#         range=[x0, x1],
#         showgrid=False),
#     yaxis2=dict(
#         domain=[0.05, 0.45],  
#         anchor="x2",
#         range=[y0, y1],
#         showgrid=False),
#     legend=dict(
#         xanchor='right',
#         yanchor='top'))

# # Copy every trace from fig1 into fig inset
# for tr in fig1.data:
#     fig.add_trace(tr, row=None, col=None)
#     fig.data[-1].update(xaxis="x2", yaxis="y2")

# Draw a rectangular border around the inset plot
# Uses "paper" coordinates so the border aligns with the inset axes domain
# fig.add_shape(
#     type="rect",
#     xref="paper", yref="paper",
#     x0=0.58, x1=0.98,
#     y0=0.05, y1=0.45,
#     line=dict(color="black", width=2.5),
#     fillcolor="rgba(0,0,0,0)",
#     opacity=1
# )

# # Draw a rectangular border around the MAIN plot area left bottom corner
# fig.add_shape(
#     type="rect",
#     x0=0, x1=60,
#     y0=0, y1=60,
#     xref="x", yref="y",
#     line=dict(color="black", width=2.5, dash="solid"),
#     fillcolor="rgba(0,0,0,0)",
#     opacity=1
# )

# # top line leading to inset
# fig.add_shape(
#     type="line",
#     x0=60, y0=60,
#     x1=0.58 * 800, y1=0.45 * 800,
#     xref="x", yref="y",
#     line=dict(color="black", width=2.5),
#     opacity=1)

# # bottom line leading to inset
# fig.add_shape(
#     type="line",
#     x0=60, y0=y0,
#     x1=0.58 * 800, y1=0.05 * 800,
#     xref="x", yref="y",
#     line=dict(color="black", width=2.5),
#     opacity=1)


fig.show()


# %%
