import pandas as pd
import numpy as np
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
pio.renderers.default = "browser"

# Create mechanism survey statistics dictionary
survey = {
    "Mechanism": ["JPMv0.2", "JPM1.1", "SmallStrato", "Superfast", "Form85", "POLLU", "GC-Hg", "E3SM", "CBM-Z", "SAPRC99", "CB05", "RACM", "RADM2", "MOZART-4", "MOZART-T1", "AMORE", "PACT-1D", "RCIM", "CIM", "CIM-var", "JAM", "GCv12.0", "GCv12.7", "GCv12.8", "GCv12.9", "GCv13.3", "GCv13.4", "GCv14.6", "CRI", "MECCA", "MCM", "Toluene", "Logan81", "CRACMM2", "CRACMM3"],
    "#Species": [11, 16, 5, 15, 12, 20, 32, 47, 67, 74, 82, 82, 59, 81, 155, 133, 167,145, 386, 394, 245, 235, 243, 258, 262, 287, 287, 353, 442, 733, 5832, 12, 38, 196, 226],
    "#Reactions": [10,13, 10, 32, 25, 25, 94, 104, 142, 211, 205, 250, 156, 196, 360, 330, 513, 379, 886, 886, 702, 725, 750, 825, 850, 903, 913, 1058, 1258, 2323, 13140, 10, 57, 531, 614],
    "Stoichiometric Invariants": [2, 5, 1, 0, 1, 3, 2, 1, 7, 1, 7, 1, 2, 0, 0, 0, 9, 0, 0, 2, 4, 6, 6, 6, 8, 10, 10, 9, 2, 9, 1, 3, 10, 4, 5],
    "Coproduction Index": [1, 1, 1, 6, 6, 4, 65, 19, 10, 13, 14, 13, 9, 17, 41, 42, 49, 85, 650, 474, 53, 136, 143, 190, 188, 185, 181, 168, 216, 456, 3558, 1, 11, 22, 35],
    "Kinetic Invariants": [1, 1, 0, 1, 0, 0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 150, 17, 0, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 1, 1, 0, 0],
    "Broken Null Cycles": [0, 0, 1, 5, 6, 4, 61, 18, 10, 13, 14, 13, 9, 17, 41, 42, 48, 85, 500, 457, 53, 135, 142, 189, 186, 183, 181, 168, 216, 456, 3558, 0, 10, 22, 35]
}

# Convert dictionary to Pandas dataframe
mech_survey = pd.DataFrame(survey)

# Create a new column for Effectuve Reactions, calculated by 
# number of total reactions "R" minus the coproduction index, or number of reactions lost to coproduction "gamma"
mech_survey["R-gamma"] = mech_survey.loc[:, "#Reactions"] - mech_survey.loc[:, "Coproduction Index"]

# Remove MCM, MECCA, CRI, GCv14.6 due to scale, remove Toluene  
mech_survey_inscale = mech_survey[~mech_survey["Mechanism"].isin(["MCM", "MECCA", "CRI", "GCv14.6", "Toluene"])]

# Filter to include only mechanisms with kinetic invariants
mech_survey_KI = mech_survey_inscale[mech_survey_inscale["Kinetic Invariants"] != 0]

# Create empty Plotly Graph Object figure
fig = go.Figure()

# Isolate specific labels, positions, and sizes per mechanism for visibility purposes
label_mechs = {"GCv13.4", "JAM", "CRACMM3", "CRACMM2", "MOZART-T1", "RCIM", "AMORE", "RACM", "SAPRC99", "CB05", "MOZART-4", "RADM2", "CBM-Z"}  
labels = [name if name in label_mechs else ""
    for name in mech_survey_inscale["Mechanism"]]
positions = ["middle right" if name in {"RCIM", "CB05", "CBM-Z"} else "middle left"
    for name in mech_survey_inscale["Mechanism"]]
sizes = [8 if name in {"Form85", "POLLU", "SmallStrato"}
         else 12 for name in mech_survey_inscale["Mechanism"]]

# Add to main scatter plot: all mechanisms, regular gray points regardless of kinetic invariants
# X-axis: # of species
# Y-axis: Effective Reactions (R - gamma)
fig.add_trace(
    go.Scatter(
        x=mech_survey_inscale["#Species"],
        y=mech_survey_inscale["R-gamma"],
        mode="markers+text",
        marker=dict(size=sizes, color="gray", opacity=1),
        text=labels,
        textposition=positions,
        textfont=dict(size=10, color='darkgray'),
        name="Mechanisms w/o Kinetic Invariants"))

# Isolate specific labels, positions, and sizes per mechanism for visibility purposes
label_mechs2 = {"GC-Hg", "Superfast", "JPM1.1", "JPMv0.2", "Logan81"} 
labels2 = ["" if name in label_mechs2 else name
    for name in mech_survey_KI["Mechanism"]]
positions1 = ["top center" if name in {"E3SM"} 
              else "bottom center" if name in {"CIM"}
              else "middle right" for name in mech_survey_KI["Mechanism"]]
sizes1 = [13 if name in {"Logan81", "Superfast", "GC-Hg", "JPM1.1", "JPMv0.2"}
          else 25 for name in mech_survey_KI["Mechanism"]]

# Add to main scatter plot: mechanisms with kinetic invariants only, purple star points
# X-axis: # of species
# Y-axis: Effective Reactions (R - gamma)
fig.add_trace(
    go.Scatter(
        x=mech_survey_KI["#Species"],
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
        name="Mechanisms w/ Kinetic Invariants"))

# Add a 1-to-1 red dotted line on the main plot
x_line = np.linspace(0, mech_survey["#Species"].max(), 100)
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
    xaxis_title="# Species",
    yaxis_title="Effective Reactions (R - γ)",
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

# Filter only the mechanisms in the crowded bottom left corner
mech_survey_crowded_KI = mech_survey[mech_survey["Mechanism"].isin(["GC-Hg", "Superfast", "JPM1.1", "JPMv0.2", "Logan81"])]
mech_survey_crowded_noKI = mech_survey[mech_survey["Mechanism"].isin(["SmallStrato", "Form85", "POLLU"])]

# Create a new empty Plotly Graph Object figure for inset plot
fig1 = go.Figure()

# Define specific text positions for visibility purposes, inset plot
positions2 = [
    "middle right" if name in {"GC-Hg", "JPM1.1"} 
    else "top center" if name in {"Superfast"}
    else "middle left" if name in {"Logan81"}
    else "bottom center" 
    for name in mech_survey_crowded_KI["Mechanism"]]

# Add to inset scatter plot: mechanisms with kinetic invariants only, purple star points
fig1.add_trace(
    go.Scatter(
        x=mech_survey_crowded_KI["#Species"],
        y=mech_survey_crowded_KI["R-gamma"],
        mode="markers+text",
        marker=dict(
            symbol="star",
            size=25,
            color="purple",
            line=dict(color="black", width=1.2)),
        text=mech_survey_crowded_KI["Mechanism"],
        textposition=positions2,
        textfont=dict(size=16),
        name="Mechanisms w/ Kinetic Invariants",
        showlegend=False))

# Define specific text positions for visibility purposes, inset plot
positions3 = [
    "middle right" if name in {"POLLU"} 
    else "middle left" if name in {"Form85"}
    else "top center" 
    for name in mech_survey_crowded_noKI["Mechanism"]]

# Add to inset scatter plot: mechanisms without kinetic invariants, regular gray points 
fig1.add_trace(
    go.Scatter(
        x=mech_survey_crowded_noKI["#Species"],
        y=mech_survey_crowded_noKI["R-gamma"],
        mode="markers+text",
        marker=dict(size=12, color="gray"),
        text=mech_survey_crowded_noKI["Mechanism"],
        textposition=positions3,
        textfont=dict(size=8, color='darkgray'),
        name="Mechanisms w/o Kinetic Invariants",
        showlegend=False))

# Add a 1-to-1 red dotted line on the inset plot
x_line = np.linspace(0, mech_survey["#Species"].max(), 100)
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

# Update additional layout properties in the inset plot
fig1.update_layout(
    xaxis=dict(tickfont=dict(size=16)),
    yaxis=dict(tickfont=dict(size=16)),
    template="simple_white",
    width=800,
    height=800,
    showlegend=False)

# Define x,y ranges for the inset plot
x0, x1 = 0, 50
y0, y1 = 0, 50

# Define the positions of the inset display area (xaxis2, yaxis2) within the main plot
fig.update_layout(
    xaxis2=dict(
        domain=[0.58, 0.98],   
        anchor="y2",
        range=[x0, x1],
        showgrid=False),
    yaxis2=dict(
        domain=[0.05, 0.45],  
        anchor="x2",
        range=[y0, y1],
        showgrid=False),
    legend=dict(
        xanchor='right',
        yanchor='top'),
    width=800, height=800)
fig.update_layout(autosize=False)

# Copy every trace from fig1 (inset plot) into fig (main plot, inset display area)
for tr in fig1.data:
    fig.add_trace(tr, row=None, col=None)
    fig.data[-1].update(xaxis="x2", yaxis="y2")

# Draw a rectangular border around the inset display area
# Uses "paper" coordinates so the border aligns with the inset axes domain
fig.add_shape(
    type="rect",
    xref="paper", yref="paper",
    x0=0.58, x1=0.98,
    y0=0.05, y1=0.45,
    line=dict(color="black", width=2.5),
    fillcolor="rgba(0,0,0,0)",
    opacity=1
)

# Draw a rectangular border around the inset plot in the main plot area, bottom left corner
fig.add_shape(
    type="rect",
    x0=0, x1=60,
    y0=0, y1=60,
    xref="x", yref="y",
    line=dict(color="black", width=2.5, dash="solid"),
    fillcolor="rgba(0,0,0,0)",
    opacity=1
)

# Add a connecting top line from the inset plot leading to the inset display area
fig.add_shape(
    type="line",
    x0=60, y0=60,
    x1=0.58 * 800, y1=0.45 * 800,
    xref="x", yref="y",
    line=dict(color="black", width=2.5),
    opacity=1)

# Add a connecting bottom line from the inset plot leading to the inset display area
fig.add_shape(
    type="line",
    x0=60, y0=y0,
    x1=0.58 * 800, y1=0.05 * 800,
    xref="x", yref="y",
    line=dict(color="black", width=2.5),
    opacity=1)

# Save figure as a PDF
pio.write_image(fig, 'survey_mechanisms.pdf', width=800, height=800) 
# Show figure in web browser
fig.show()
