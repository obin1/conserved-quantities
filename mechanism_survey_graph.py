import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes

survey = {
    "Mechanism": ["JPMv0.2", "JPM1.1", "SmallStrato", "Superfast", "Form 1985", "POLLU", "GC-Hg", "E3SM", "CBM-Z", "SAPRC99", "CB05", "RACM", "RADM2", "MOZART", "MOZART-T1", "AMORE", "PACT-1D", "RICM", "CIM", "CIM-var", "JAM", "GCv12.0", "GCv12.7", "GCv12.8", "GCv12.9", "GCv13.3", "GCv13.4", "GCv14.6", "CRI", "MECCA", "MCM", "Daniel"],
    "#Species": [11, 16, 5, 15, 12, 20, 32, 47, 67, 74, 82, 82, 59, 81, 155, 133, 167,145, 386, 394, 245, 235, 243, 258, 262, 287, 287, 353, 442, 733, 5832, 12],
    "#Reactions": [10,13, 10, 32, 25, 25, 94, 104, 142, 211, 205, 250, 156, 196, 360, 330, 513, 379, 886, 886, 702, 725, 750, 825, 850, 903, 913, 1058, 1258, 2323, 13140, 10],
    "Stoichiometric Invariants": [2, 5, 1, 0, 1, 3, 2, 1, 7, 1, 7, 1, 2, 0, 0, 0, 9, 0, 0, 2, 4, 6, 6, 6, 8, 10, 10, 9, 2, 9, 1, 3],
    "Coproduction Index": [1, 1, 1, 6, 6, 4, 65, 19, 10, 13, 14, 13, 9, 17, 41, 42, 49, 85, 650, 474, 53, 136, 143, 190, 188, 185, 181, 168, 216, 456, 3558, 1],
    "Kinetic Invariants": [1, 1, 0, 1, 0, 0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 150, 17, 0, 1, 1, 1, 2, 2, 0, 0, 0, 0, 0, 1],
    "Broken Null Cycles": [0, 0, 1, 5, 6, 4, 61, 18, 10, 13, 14, 13, 9, 17, 41, 42, 48, 85, 500, 457, 53, 135, 142, 189, 186, 183, 181, 168, 216, 456, 3558, 0]
}

mech_survey = pd.DataFrame(survey)

mech_survey["R-gamma"] = mech_survey.loc[:, "#Reactions"] - mech_survey.loc[:, "Coproduction Index"]

mech_survey = mech_survey[mech_survey["Mechanism"] != "MCM"]
mech_survey = mech_survey[mech_survey["Mechanism"] != "MECCA"]
mech_survey = mech_survey[mech_survey["Mechanism"] != "Daniel"]

mech_survey2 = mech_survey[mech_survey["Kinetic Invariants"] != 0]


plt.scatter(mech_survey.loc[:, "#Species"], mech_survey.loc[:, "R-gamma"])
plt.scatter(mech_survey2.loc[:, "#Species"], mech_survey2.loc[:, "R-gamma"], marker="*", s=150, linewidths=1.2, edgecolors="black", color='green')

x = np.linspace(0, mech_survey[["#Species"]].max()[0], 100) 
y = x
plt.plot(x, y, color='red', linestyle='--')
plt.fill_between(x, 0, y, alpha=0.2, color='purple')

#%%

fig, ax = plt.subplots(figsize=(10, 8))

# Scatter: mech_survey
ax.scatter(
    mech_survey["#Species"],
    mech_survey["R-gamma"],
    s=64,                     # Plotly size=8 ≈ matplotlib s=64
    label="Mechanisms w/o Kinetic Invariants"
)

# Scatter: mech_survey2 (stars + labels)
ax.scatter(
    mech_survey2["#Species"],
    mech_survey2["R-gamma"],
    marker="*",
    s=250,                    # Plotly size=16 ≈ matplotlib s=250
    color="green",
    edgecolors="black",
    linewidths=1.2,
    label="Mechanisms w/ Kinetic Invariants"
)

# Text labels for starred points
for _, row in mech_survey2.iterrows():
    ax.annotate(
        row["Mechanism"],
        (row["#Species"], row["R-gamma"]),
        textcoords="offset points",
        xytext=(15, -5),
        ha="left",
        fontsize=10
    )

# 1–1 line
x = np.linspace(0, mech_survey["#Species"].max(), 100)
ax.plot(x,x,linestyle="--",color="red",label="1–1 line")

# Shaded region below 1–1 line
ax.fill_between(x,0,x,color="purple",alpha=0.2)

# Labels & styling
ax.set_xlabel("# Species")
ax.set_ylabel("R − γ")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %%
fig = plt.figure(figsize=(12, 10))

bax = brokenaxes(ylims=((0, 350), (600, 750)), hspace=0.05, fig=fig)

# Scatter: mech_survey
bax.scatter(
    mech_survey["#Species"],
    mech_survey["R-gamma"],
    label="Mechanisms w/o Kinetic Invariants",
    s=50
)

# Scatter: mech_survey2 (stars + labels)
bax.scatter(
    mech_survey2["#Species"],
    mech_survey2["R-gamma"],
    marker="*",
    s=150,
    linewidths=1.2,
    edgecolors="black",
    color="green",
    label="Mechanisms w/ Kinetic Invariants"
)

# Annotate stars with mechanism names
for _, row in mech_survey2.iterrows():
    bax.annotate(
        row["Mechanism"],
        (row["#Species"], row["R-gamma"]),
        textcoords="offset points",
        xytext=(15, -5),
        fontsize=12)

# 1–1 line
x = np.linspace(0, mech_survey["#Species"].max()+50, 100)
bax.plot(x, x, color="red", linestyle="--", label="1–1 line")
x2 = np.linspace(350, 500, 100)
y2 = np.linspace(590, 750, 100)
bax.plot(x2, y2, color="red", linestyle="--")

# Purple shading below the line
bax.fill_between(x, 0, x, alpha=0.2, color="purple")
bax.fill_between(x2, 600, y2, alpha=0.2, color="purple")

# Labels and legend
bax.set_xlabel("# Species")
bax.set_ylabel("R − γ")
bax.legend(loc="upper left")

plt.show()
# %%
fig = plt.figure(figsize=(12, 10))

bax = brokenaxes(ylims=((0, 350), (600, 1100)), hspace=0.05, fig=fig)

# Scatter: mech_survey
bax.scatter(
    mech_survey["#Species"],
    mech_survey["R-gamma"],
    label="Mechanisms w/o Kinetic Invariants",
    s=50
)

# Scatter: mech_survey2 (stars + labels)
bax.scatter(
    mech_survey2["#Species"],
    mech_survey2["R-gamma"],
    marker="*",
    s=150,
    linewidths=1.2,
    edgecolors="black",
    color="green",
    label="Mechanisms w/ Kinetic Invariants"
)

# Annotate stars with mechanism names
for _, row in mech_survey2.iterrows():
    bax.annotate(
        row["Mechanism"],
        (row["#Species"], row["R-gamma"]),
        textcoords="offset points",
        xytext=(15, -5),
        fontsize=12)

# 1–1 line
x = np.linspace(0, 600, 100)
bax.plot(x, x, color="red", linestyle="--", label="1–1 line")
x2 = np.linspace(350, 760, 100)
y2 = np.linspace(590, 1000, 100)
bax.plot(x2, y2, color="red", linestyle="--")

# Purple shading below the line
bax.fill_between(x, 0, x, alpha=0.2, color="purple")
bax.fill_between(x2, 600, y2, alpha=0.2, color="purple")

# Labels and legend
bax.set_xlabel("# Species")
bax.set_ylabel("R − γ")
bax.legend(loc="upper left")

plt.show()
# %%
import plotly.graph_objects as go
import numpy as np
import plotly.io as pio
pio.renderers.default = "browser"

#%%

# Create figure
fig = go.Figure()

# Scatter: mech_survey (regular points)
fig.add_trace(
    go.Scatter(
        x=mech_survey["#Species"],
        y=mech_survey["R-gamma"],
        mode="markers",
        marker=dict(size=8, color="blue"),
        name="Mechanisms w/o Kinetic Invariants"
    )
)

# Scatter: mech_survey2 (stars)
fig.add_trace(
    go.Scatter(
        x=mech_survey2["#Species"],
        y=mech_survey2["R-gamma"],
        mode="markers+text",
        marker=dict(
            symbol="star",
            size=16,
            color="green",
            line=dict(color="black", width=1.2)
        ),
        text=mech_survey2["Mechanism"],
        textposition="middle right",
        name="Mechanisms w/ Kinetic Invariants"
    )
)

# 1-1 line
x_line = np.linspace(0, mech_survey["#Species"].max(), 100)
y_line = x_line
fig.add_trace(
    go.Scatter(
        x=x_line,
        y=y_line,
        mode="lines",
        line=dict(color="red", dash="dash"),
        name="1–1 line"
    )
)

# Fill under 1-1 line
fig.add_trace(
    go.Scatter(
        x=np.concatenate([x_line, x_line[::-1]]),
        y=np.concatenate([np.zeros_like(x_line), y_line[::-1]]),
        fill="toself",
        fillcolor="rgba(128,0,128,0.2)",
        line=dict(color="rgba(255,255,255,0)"),
        hoverinfo="skip",
        showlegend=False
    )
)

# Layout
fig.update_layout(
    xaxis_title="# Species",
    yaxis_title="R − γ",
    template="simple_white",
    width=1000,
    height=800
)

fig.show()
# %%
