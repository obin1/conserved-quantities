import pandas as pd
import numpy as np
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
from sklearn.linear_model import LinearRegression
pio.renderers.default = "browser"

# Read CSV file of mechanism survey
mech_survey = pd.read_csv("mechanism_survey.csv").dropna()

# Create a new column for Effectuve Reactions, calculated by 
# number of total reactions "R" minus the coproduction index, or number of reactions lost to coproduction "gamma"
mech_survey["Reactions"] = pd.to_numeric(mech_survey["Reactions"], errors="coerce")
mech_survey["R-gamma"] = mech_survey.loc[:, "Reactions"] - mech_survey.loc[:, "Coproduction Index"]

# Remove MCM, MECCA, CRI, GCv14.6 due to scale, remove Toluene  
mech_survey_inscale = mech_survey[~mech_survey["Short Name"].isin(["MCM", "MECCA", "CRI", "GCv14.6", "MCM-PRAM"])]

# Filter to include only mechanisms with kinetic invariants
mech_survey_KI = mech_survey_inscale[mech_survey_inscale["Kinetic Invariants"] != 0]

# Create empty Plotly Graph Object figure
fig = go.Figure()

# Isolate specific labels, positions, and sizes per mechanism for visibility purposes
label_mechs = {"GCv13.4", "JAM", "CRACMM3", "CRACMM2", "MOZART-T1", "RCIM", "AMORE", "RACM", "SAPRC99", "CB05", "MOZART-4", "RADM2", "CBM-Z"}  
labels = [name if name in label_mechs else ""
    for name in mech_survey_inscale["Short Name"]]
positions = ["middle right" if name in {"RCIM", "CB05", "CBM-Z"} else "middle left"
    for name in mech_survey_inscale["Short Name"]]
sizes = [8 if name in {"Form85", "POLLU", "SmallStrato"}
         else 12 for name in mech_survey_inscale["Short Name"]]

# Add to main scatter plot: all mechanisms, regular gray points regardless of kinetic invariants
# X-axis: # of species
# Y-axis: Effective Reactions (R - gamma)
fig.add_trace(
    go.Scatter(
        x=mech_survey_inscale["Species"],
        y=mech_survey_inscale["Reactions"],
        mode="markers+text",
        marker=dict(size=sizes, color="gray", opacity=1),
        text=labels,
        textposition=positions,
        textfont=dict(size=10, color='darkgray'),
        name="without Kinetic Invariants"))

# Isolate specific labels, positions, and sizes per mechanism for visibility purposes
label_mechs2 = {"GC-Hg", "Superfast", "JPM1.1", "JPMv0.2", "Logan81"} 
labels2 = ["" if name in label_mechs2 else name
    for name in mech_survey_KI["Short Name"]]
positions1 = ["top center" if name in {"E3SM"} 
              else "bottom center" if name in {"CIM"}
              else "middle right" for name in mech_survey_KI["Short Name"]]
sizes1 = [13 if name in {"Logan81", "Superfast", "GC-Hg", "JPM1.1", "JPMv0.2"}
          else 25 for name in mech_survey_KI["Short Name"]]

# Add to main scatter plot: mechanisms with kinetic invariants only, purple star points
# X-axis: # of species
# Y-axis: Effective Reactions (R - gamma)
fig.add_trace(
    go.Scatter(
        x=mech_survey_KI["Species"],
        y=mech_survey_KI["Reactions"],
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

# all mechanisms
lr = LinearRegression(fit_intercept=False)
x = mech_survey[["Species"]]
y = mech_survey[["Reactions"]]
lr.fit(x, y)
m = lr.coef_[0]
x_line_reg = np.linspace(0, mech_survey["Species"].max(), 100)
y_line_reg = x_line_reg*m
fig.add_trace(
    go.Scatter(
        x=x_line_reg,
        y=y_line_reg,
        mode="lines",
        line=dict(color="blue"),
        name=f'{m[0]:.2f} line',
        showlegend=True))

# no CIM & CIM-var
mech_survey_exclude = mech_survey[~mech_survey["Short Name"].isin(["CIM", "CIM-var", "MCM", "MECCA", "CRI", "GCv14.6", "MCM-PRAM"])]

lr2 = LinearRegression(fit_intercept=False)
x2 = mech_survey_exclude[["Species"]]
y2 = mech_survey_exclude[["Reactions"]]
lr2.fit(x2, y2)
m2 = lr2.coef_[0]
x_line_reg2 = np.linspace(0, mech_survey_exclude["Species"].max(), 100)
y_line_reg2 = x_line_reg2*m2
fig.add_trace(
    go.Scatter(
        x=x_line_reg2,
        y=y_line_reg2,
        mode="lines",
        line=dict(color="orange"),
        name=f'{m2[0]:.2f} line',
        showlegend=True))

# Add a 1-to-1 red dotted line on the main plot
x_line = np.linspace(0, mech_survey["Species"].max(), 100)
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
    yaxis_title="Reactions",
    xaxis=dict(
        title_font=dict(size=20), 
        tickfont=dict(size=16)),
    yaxis=dict(
        title_font=dict(size=20),  
        tickfont=dict(size=16)    ),
    template="simple_white",
    width=800,
    height=800,
    legend=dict(
        orientation="h",
        x=0.85,
        xanchor="center",
        y=0.99,
        yanchor="top",
    ),)
fig.update_xaxes(range=[0, 800])
fig.update_yaxes(range=[0, 800])


# Save figure as a PDF
pio.write_image(fig, 'survey_mechanisms_regression.pdf', width=800, height=800) 
# Show figure in web browser
fig.show()
