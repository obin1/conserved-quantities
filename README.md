# conserved-quantities

Repository for analysis of conserved quantities in chemical mechanisms.

Contents
- Julia scripts for computing conserved cycles (e.g. `ConservedCycles.jl`, `DrivingForceInvariant.jl`).
- Example data files (biadjacency matrices, CSVs).

Quick start
1. Open this folder in VS Code or a Julia REPL.
2. Run the main scripts with Julia 1.8+:

```julia
# from project root
include("ConservedCycles.jl")
# or run your own REPL workflow
```

Notes
- This repository was initialized locally; add a remote if you want to push to GitHub/GitLab:

```bash
git remote add origin <your-remote-url>
git push -u origin main
```

License
This repository is published under the MIT License. See `LICENSE`.
