import pandas as pd
import numpy as np
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
pio.renderers.default = "browser"

# Read CSV file of mechanism survey
mech_survey = pd.read_csv("mechanism_survey.csv").dropna()

# Create a new column for Effectuve Reactions, calculated by 
# number of total reactions "R" minus the coproduction index, or number of reactions lost to coproduction "gamma"
mech_survey["Reactions"] = pd.to_numeric(mech_survey["Reactions"], errors="coerce")
mech_survey["R-gamma"] = mech_survey.loc[:, "Reactions"] - mech_survey.loc[:, "Coproduction Index"]

# Remove MCM, MECCA, CRI, GCv14.6 due to scale, remove Toluene  
# mech_survey_inscale = mech_survey[~mech_survey["Short Name"].isin(["MCM", "MECCA", "CRI", "GCv14.6"])]

# Filter to include only mechanisms with kinetic invariants
mech_survey_KI = mech_survey[mech_survey["Kinetic Invariants"] != 0]

# Create empty Plotly Graph Object figure
fig = go.Figure()

# Isolate specific labels, positions, and sizes per mechanism for visibility purposes
label_mechs = {"Form85", "POLLU", "SmallStrato", "MECCA", "MCM", "CRI"}  
labels = [name if name in label_mechs else ""
    for name in mech_survey["Short Name"]]
positions = ["middle right" if name in {"RCIM", "CB05", "CBM-Z", "MOZART-4"} 
             else "top center" if name in {"SmallStrato"}
             else "bottom center" if name in {"POLLU"}
             else "middle left" if name in {"Form85"}
             else "middle left"
    for name in mech_survey["Short Name"]]
sizes = [12 if name in label_mechs
         else 8 for name in mech_survey["Short Name"]]

# Add to main scatter plot: all mechanisms, regular gray points regardless of kinetic invariants
# X-axis: # of species
# Y-axis: Effective Reactions (R - gamma)
fig.add_trace(
    go.Scatter(
        x=mech_survey["Species"],
        y=mech_survey["R-gamma"],
        mode="markers+text",
        # mode="text",
        marker=dict(size=sizes, color="gray", opacity=1),
        text=labels,
        textposition=positions,
        textfont=dict(size=10, color='darkgray'),
        name="without Kinetic Invariants"))

# Isolate specific labels, positions, and sizes per mechanism for visibility purposes
# label_mechs2 = {"GC-Hg", "Superfast", "JPM1.1", "JPMv0.2", "Logan81"} 
label_mechs2 = {"JPMv0.2", "JPMv1.1", "E3SM", "GC-Hg", "Logan81", "Superfast", "MCM-PRAM"}
labels2 = [name if name in label_mechs2 else "" for name in mech_survey_KI["Short Name"]]
positions1 = ["bottom center" if name in {"CIM", "JPMv0.2"}
              else "middle left" if name in {"Superfast", "E3SM"}
              else "middle right" for name in mech_survey_KI["Short Name"]]
sizes1 = [25 if name in label_mechs2
          else 15 for name in mech_survey_KI["Short Name"]]

# Add to main scatter plot: mechanisms with kinetic invariants only, purple star points
# X-axis: # of species
# Y-axis: Effective Reactions (R - gamma)
fig.add_trace(
    go.Scatter(
        x=mech_survey_KI["Species"],
        y=mech_survey_KI["R-gamma"],
        mode="markers+text",
        marker=dict(
            symbol="star",
            size=sizes1,
            color="purple",
            opacity=1,
            line=dict(color="black", width=1.2)),
        text=labels2,
        textposition=positions1,
        textfont=dict(size=16),
        name="with Kinetic Invariants"))

# Add a 1-to-1 red dotted line on the main plot
# x_line = np.linspace(0, mech_survey["Species"].max(), 100)
x_line = np.logspace(
    np.log10(1),
    np.log10(100000),
    100
)
y_line = x_line
fig.add_trace(
    go.Scatter(
        x=x_line,
        y=y_line,
        mode="lines",
        line=dict(color="red", dash="dash"),
        name="1–1 line",
        showlegend=False))

# Fill the area under the 1-to-1 line in the main plot
fig.add_trace(
    go.Scatter(
        x=np.concatenate([x_line, x_line[::-1]]),
        y=np.concatenate([np.zeros_like(x_line), y_line[::-1]]),
        fill="toself",
        fillcolor="rgba(128,0,128,0.2)",
        line=dict(color="rgba(255,255,255,0)"),
        hoverinfo="skip",
        showlegend=False))

# Update additional layout properties in the main plot
fig.update_layout(
    xaxis_title="Species",
    yaxis_title="Effective Reactions",
    xaxis=dict(
        title_font=dict(size=20), 
        tickfont=dict(size=16)),
    yaxis=dict(
        title_font=dict(size=20),  
        tickfont=dict(size=16)),
    template="simple_white",
    width=2500,
    height=2000,
    margin=dict(l=90, r=0, t=0, b=0),
    legend=dict(
        orientation="h",
        x=0.85,
        xanchor="center",
        y=0.92,
        yanchor="top",
        font=dict(
            size=15 
        )
    ),
    autosize=False
)

fig.update_xaxes(
    type="log",
    tickmode="array",
    tickvals=[1, 10, 100, 1000, 10000, 100000],
    ticktext=["10⁰", "10¹", "10²", "10³", "10⁴", "10⁵"], constrain="domain",
    autorange=False, 
    range=[0, 5]
)

fig.update_yaxes(
    type="log",
    tickmode="array",
    tickvals=[1, 10, 100, 1000, 10000, 100000],
    ticktext=["10⁰", "10¹", "10²", "10³", "10⁴", "10⁵"],
    scaleanchor="x",
    scaleratio=1, constrain="domain",
    autorange=False, 
    range=[0, 5]
)


# Filter only the mechanisms in the crowded bottom left corner
mech_survey_crowded_KI = mech_survey[mech_survey["Short Name"].isin(["CIM", "CIM-var", "PACT-1D", "GCv12.0", "GCv12.7", "GCv12.8", "GCv12.9", "GCv13.3"])]
mech_survey_crowded_noKI = mech_survey[mech_survey["Short Name"].isin(["RADM2", "MOZART-4", "CB05", "RACM", "SAPRC99", "AMORE", "RCIM", "MOZART-T1", "CRACMM2", "CRACMM3", "JAM", "GCv13.4", "GCv14.6", "CRI"])]


# Create a new empty Plotly Graph Object figure for inset plot
fig1 = go.Figure()

# Define specific text positions for visibility purposes, inset plot
positions2 = [
    "bottom center" if name in {"CIM-var"} 
    else "bottom right" if name in {"GCv12.0"}
    else "middle right" 
    for name in mech_survey_crowded_KI["Short Name"]]

# Add to inset scatter plot: mechanisms with kinetic invariants only, purple star points
fig1.add_trace(
    go.Scatter(
        x=mech_survey_crowded_KI["Species"],
        y=mech_survey_crowded_KI["R-gamma"],
        mode="markers+text",
        marker=dict(
            symbol="star",
            size=25,
            color="purple",
            line=dict(color="black", width=1.2)),
        text=mech_survey_crowded_KI["Short Name"],
        textposition=positions2,
        # textposition="middle right",
        textfont=dict(size=11),
        name="with Kinetic Invariants",
        showlegend=False))

# Define specific text positions for visibility purposes, inset plot
positions3 = [
    "middle right" if name in {"RADM2", "MOZART-4", "MOZART-T1", "RCIM", "SAPRC99", "RACM"} 
    else "middle left" if name in {"CRACMM2", "CRACMM3", "CB05", "GCv14.6", "CRI", "AMORE"}
    # else "bottom right" if name in {"AMORE"}
    else "top left" if name in {"GCv13.4", "JAM"}
    # else "top right" if name in {"CB05"}
    else "top center" 
    for name in mech_survey_crowded_noKI["Short Name"]]

# Add to inset scatter plot: mechanisms without kinetic invariants, regular gray points 
fig1.add_trace(
    go.Scatter(
        x=mech_survey_crowded_noKI["Species"],
        y=mech_survey_crowded_noKI["R-gamma"],
        mode="markers+text",
        marker=dict(size=9, color="gray"),
        text=mech_survey_crowded_noKI["Short Name"],
        textposition=positions3,
        textfont=dict(size=8, color='darkgray'),
        name="without Kinetic Invariants",
        showlegend=False))

# Add a 1-to-1 red dotted line on the inset plot
x_line = np.logspace(
    np.log10(1),
    np.log10(100000),
    100
)
y_line = x_line

fig1.add_trace(
    go.Scatter(
        x=x_line,
        y=y_line,
        mode="lines",
        line=dict(color="red", dash="dash"),
        name="1–1 line",
        showlegend=False))

# Fill the area under the 1-to-1 line in the inset plot
fig1.add_trace(
    go.Scatter(
        x=np.concatenate([x_line, x_line[::-1]]),
        y=np.concatenate([np.zeros_like(x_line), y_line[::-1]]),
        fill="toself",
        fillcolor="rgba(128,0,128,0.2)",
        line=dict(color="rgba(255,255,255,0)"),
        hoverinfo="skip",
        showlegend=False))

# Define x,y ranges for the inset plot
x0, x1 = 1, 550
y0, y1 = 125, 950

# Update additional layout properties in the inset plot
# Define the positions of the inset display area (xaxis2, yaxis2) within the main plot
fig.update_layout(
    xaxis2=dict(
        domain=[0.65, 0.98],   
        anchor="y2",
        range=[x0, x1],
        showgrid=False,
        ),
    yaxis2=dict(
        domain=[0.12, 0.55],  
        anchor="x2",
        range=[y0, y1],
        showgrid=False,
        ),
    legend=dict(
        xanchor='right',
        yanchor='top'),
    width=800, height=800)

# Copy every trace from fig1 (inset plot) into fig (main plot, inset display area)
for tr in fig1.data:
    fig.add_trace(tr, row=None, col=None)
    fig.data[-1].update(xaxis="x2", yaxis="y2")

# Draw a rectangular border around the inset display area
# Uses "paper" coordinates so the border aligns with the inset axes domain
fig.add_shape(
    type="rect",
    xref="paper", yref="paper",
    x0=0.65, x1=0.98,
    y0=0.12, y1=0.55,
    line=dict(color="black", width=1.5),
    fillcolor="rgba(0,0,0,0)",
    opacity=1
)

# Draw a rectangular border around the inset plot in the main plot area
fig.add_shape(
    type="rect",
    x0=40, x1=550,
    y0=125, y1=950,
    xref="x", yref="y",
    line=dict(color="black", width=1.5, dash="solid"),
    fillcolor="rgba(0,0,0,0)",
    opacity=1
)

# Add a connecting top line from the inset plot leading to the inset display area
# fig.add_shape(
#     type="line",
#     x0=300, y0=950,
#     x1=1000, y1=1500,
#     xref="x", yref="y",
#     line=dict(color="black", width=1.5),
#     opacity=1)

# # Add a connecting bottom line from the inset plot leading to the inset display area
# fig.add_shape(
#     type="line",
#     x0=300, y0=120,
#     x1=1000, y1=1,
#     xref="x", yref="y",
#     line=dict(color="black", width=1.5),
#     opacity=1)

# Save figure as a PDF
pio.write_image(fig, 'survey_mechanisms_log.pdf', width=800, height=800) 
# Show figure in web browser
fig.show()
