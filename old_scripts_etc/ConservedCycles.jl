using LinearAlgebra
using RowEchelon
using DelimitedFiles
using SparseArrays

function conservationQR(N)
#=
This function takes in a biadjacency matrix (stoichiometric matrix)
of a species-reaction network, where rows correspond to species and
columns correspond to rows.
e.g. A in equation 2 of Sturm and Wexler (2020) or N in equation 1
of Vallabhajosyula et al. (2006).
The output of this function is conservation laws found using the QR 
factorization method from Vallabhajosyula et al. (2006). The matrix
Γ specifies the conserved quantities.

To reproduce the example from Vallabhajosyula et al. (2006)
N = [-1  0  1;
      1  0 -1;
      1 -1  0;
      0  1 -1 ]
Note that R will be permuted from Vallabhajosyula et al. (2006)
The permutations p will also be different
This is fine: the same conserved cycles hold
=#
    # Perform Householder QR factorization on the transpose of N
    NT = transpose(N)[:,:]
    Q,R, p = qr(NT, ColumnNorm())

    # Put R matrix in reduced row echelon form and round to zero
    Rstar = rref(R)
    Rstar[abs.(Rstar).<1e-10] .= 0

    # L0 is the transpose of M, which is the nonzero, non-identity block
    L0 = transpose(Rstar[1:rank(Rstar),rank(Rstar)+1:end])

    # Γ is L0, with identity horizontally concatenated
    Γ = hcat(-L0,I(size(L0)[1]))*I(size(N)[1])[p,:]

    return transpose(Γ)
end

function print_conserved_cycles(spc_names, Γ; tol=1e-3)
#= 
This function takes in a conservation matrix and a species list 
and outputs the conserved cycles as a list of strings.    
=#
    conserved_cycles = []
    

    for j = 1:size(Γ)[2]
        Γ[:,j] = map(x -> abs(x) < tol ? 0 : x, Γ[:,j])
        spc_list = spc_names[Γ[:,j] .!= 0]
        spc_coeffs = string.(round.(Γ[Γ[:,j] .!= 0, j], digits=5))
        spc_coeffs[spc_coeffs .== "1.0"] .= " "
        push!(conserved_cycles, join(spc_coeffs .* spc_list, " + "))
    end
    return conserved_cycles
end

# Reproducing the example from Vallabhajosyula et al. (2006)
N = [-1  0  1;
      1  0 -1;
      1 -1  0;
      0  1 -1 ]

Γtest = conservationQR(N)
# Note that R will be permuted from https://doi.org/10.1093/bioinformatics/bti800
# The permutations p will also be different
# This is fine: the same conserved cycles hold

# Application to the Leighton relationship
# R1: NO2 + hv → NO + O
# R2: O   + O2 → O3
# R3: O3  + NO → NO2
# The first 3 species are ordered like JPM: 
# O3, NO, NO2, then O, then O2
A_leight =   [0 -1  1;
              1  0 -1;
             -1  0  1;
              1 -1  0;
              0 -1  0]

Γ_leight = conservationQR(A_leight)




# Application to the Julia photochemical model
B =           [0  1 -1  0  0  0  0  0  0  0;
               1  0 -1  0  0  0 -1  0  0  0;
              -1  0  1  0  0  0  1 -1  0  0;
               0  0  0 -1 -1 -1  0  0  0  0;
               0  0  0  2  0  1 -1  0  0  1;
               0  0  0  0  0  0  0  0 -1 -1;
               0  0  0  0  0 -1  1 -1  2 -1;
               1 -1  0  0  0  0  0  0  0  0;
               0  0  0  0  0  0  0  1  0  0;
               0  0  0  1  1  1  0  0  0  0;
               0  0  0  0  1  0  0  0  0  0]
Γ = conservationQR(B)
jpm_spc_names = ["O3","NO","NO2","HCHO","HO2.","HO2H","OH.","O","HNO3","CO","H2"]
jpm_conserved_cycles = print_conserved_cycles(jpm_spc_names,Γ)

# Application to Super Fast Chemistry
# Note: no conserved cycles!
A_sup = readdlm("/Users/psturm/Desktop/KPP-playground/superfast/superfast_BiadjacencyMatrix.csv",',',skipstart=26)
sup_spc_names = ["" for i = 1:maximum(A_sup[:,2])]
sup_spc_names[A_sup[:,2]] .= A_sup[:,1]
A_sup = Float64.(A_sup[:,2:end])
A_sup = sparse(A_sup[:,1],A_sup[:,2],A_sup[:,3])
Γ_sup = conservationQR(A_sup)

# Application to Small Strato
# A_sml = readdlm("/Users/psturm/Desktop/KPP-playground/small_strato_example/small_strato_BiadjacencyMatrix.csv",',',skipstart=26)
# sml_spc_names = ["" for i = 1:maximum(A_sml[:,2])]
# sml_spc_names[A_sml[:,2]] .= A_sml[:,1]
# A_sml = Float64.(A_sml[:,2:end])
# A_sml = sparse(A_sml[:,1],A_sml[:,2],A_sml[:,3])
# Γ_sml = conservationQR(A_sml)
# print_conserved_cycles(sml_spc_names,Γ_sml)

# Application to RACM
# Found here: https://mbees.med.uni-augsburg.de/gitlab/mbees/boxmox/-/blob/master/examples/RACM.kpp
A_racm = readdlm("/Users/psturm/Desktop/KPP-playground/RACM/RACM_BiadjacencyMatrix.csv",',',skipstart=26)
racm_spc_names = ["" for i = 1:maximum(A_racm[:,2])]
racm_spc_names[A_racm[:,2]] .= A_racm[:,1]
A_racm = Float64.(A_racm[:,2:end])
A_racm = sparse(A_racm[:,1],A_racm[:,2],A_racm[:,3])
Γ_racm = conservationQR(A_racm)
print_conserved_cycles(racm_spc_names,Γ_racm)

# Application to CBMZ
# Found here: https://mbees.med.uni-augsburg.de/gitlab/mbees/boxmox/-/tree/master/examples
A_cbmz = readdlm("/Users/psturm/Desktop/KPP-playground/CBMZ/CBMZ_BiadjacencyMatrix.csv",',',skipstart=26)
cbmz_spc_names = ["" for i = 1:maximum(A_cbmz[:,2])]
cbmz_spc_names[A_cbmz[:,2]] .= A_cbmz[:,1]
A_cbmz = Float64.(A_cbmz[:,2:end])
A_cbmz = sparse(A_cbmz[:,1],A_cbmz[:,2],A_cbmz[:,3])
Γ_cbmz = conservationQR(A_cbmz)
# divide each column by its smallest nonzero element
for j = 1:size(Γ_cbmz)[2]
    min_val = minimum(abs.(Γ_cbmz[:,j][Γ_cbmz[:,j] .!= 0]))
    Γ_cbmz[:,j] = Γ_cbmz[:,j] ./ min_val
end
print_conserved_cycles(cbmz_spc_names,Γ_cbmz)
# how many nonzero elements are there in each column?
sum(abs.(Γ_cbmz) .> 1e-10,dims=1)
# OKAY ALSO there are 10 missing indices in the species list, so these species get no spc_name
cbmz =readdlm("/Users/psturm/Desktop/KPP-playground/CBMZ/CBMZ_BiadjacencyMatrix.csv",',',skipstart=26)
unique(cbmz[:,2])

# Application to MOZART
# Found here: http://meteo.edu.vn/~trungnq/Download/Save_WRFCHEM_for_Sri/WRFV3/chem/KPP/mechanisms/mozart/
A_moz = readdlm("/Users/psturm/Desktop/KPP-playground/mozart/mozart_BiadjacencyMatrix.csv",',',skipstart=26)
moz_spc_names = ["" for i = 1:maximum(A_moz[:,2])]
moz_spc_names[A_moz[:,2]] .= A_moz[:,1]
A_moz = Float64.(A_moz[:,2:end])
A_moz = sparse(A_moz[:,1],A_moz[:,2],A_moz[:,3])
Γ_moz = conservationQR(A_moz)

# Application to MOZART-T1
# Found here: https://mbees.med.uni-augsburg.de/gitlab/mbees/boxmox/-/blob/master/examples/MOZART_T1.kpp
A_mt1 = readdlm("/Users/psturm/Desktop/KPP-playground/mozart_t1/MOZART_T1_BiadjacencyMatrix.csv",',',skipstart=26)
mt1_spc_names = ["" for i = 1:maximum(A_mt1[:,2])]
mt1_spc_names[A_mt1[:,2]] .= A_mt1[:,1]
A_mt1 = Float64.(A_mt1[:,2:end])
A_mt1 = sparse(A_mt1[:,1],A_mt1[:,2],A_mt1[:,3])
Γ_mt1 = conservationQR(A_mt1)

# Application to SAPRC99
# From WRF github: https://github.com/wrf-model/WRF/tree/master/chem/KPP/mechanisms/saprc99
# A_sap = readdlm("/Users/psturm/Desktop/KPP-playground/saprc99/saprc99_BiadjacencyMatrix.csv",',',skipstart=26)
# sap_spc_names = ["" for i = 1:maximum(A_sap[:,2])]
# sap_spc_names[A_sap[:,2]] .= A_sap[:,1]
# A_sap = Float64.(A_sap[:,2:end])
# A_sap = sparse(A_sap[:,1],A_sap[:,2],A_sap[:,3])
# Γ_sap = conservationQR(A_sap)
# print_conserved_cycles(sap_spc_names,Γ_sap)

# Application to CB05-TUCl, EPA
# Found here: https://mbees.med.uni-augsburg.de/gitlab/mbees/boxmox/-/tree/master/examples
A_cb05 = readdlm("/Users/psturm/Desktop/KPP-playground/cb05_tucl_epa/CB05TUCl_EPA_BiadjacencyMatrix.csv",',',skipstart=26)
cb05_spc_names = ["" for i = 1:maximum(A_cb05[:,2])]
cb05_spc_names[A_cb05[:,2]] .= A_cb05[:,1]
A_cb05 = Float64.(A_cb05[:,2:end])
A_cb05 = sparse(A_cb05[:,1],A_cb05[:,2],A_cb05[:,3])
Γ_cb05 = conservationQR(A_cb05)
print_conserved_cycles(cb05_spc_names,Γ_cb05)

# Application to Caltech Isoprene Mechanism, Reduced Plus
# From here: https://data.caltech.edu/records/x88rk-wca37
# Reduced Plus ran after renaming CH2O to HCHO in 3 places
A_iso = readdlm("/Users/psturm/Desktop/KPP-playground/ReducedPlus/isoprene_reduced_plus_v5_BiadjacencyMatrix.csv",',',skipstart=26)
iso_spc_names = ["" for i = 1:maximum(A_iso[:,2])]
iso_spc_names[A_iso[:,2]] .= A_iso[:,1]
A_iso = Float64.(A_iso[:,2:end])
A_iso = sparse(A_iso[:,1],A_iso[:,2],A_iso[:,3])
Γ_iso = conservationQR(A_iso)
print_conserved_cycles(iso_spc_names,Γ_iso)

# Application to GEOS-Chem v14.0
A_gc = readdlm("/Users/psturm/Desktop/Twilight_KPP/CheckKPPStandalone/gckpp_BiadjacencyMatrix.csv",',',skipstart=26)
gc_spc_names = ["" for i = 1:maximum(A_gc[:,2])]
gc_spc_names[A_gc[:,2]] .= A_gc[:,1]
A_gc = Float64.(A_gc[:,2:end])
A_gc = sparse(A_gc[:,1],A_gc[:,2],A_gc[:,3])
Γ_gc = conservationQR(A_gc)
print_conserved_cycles(gc_spc_names,Γ_gc)
L = nullspace(Matrix(A_gc)')

# Is NOx family conservation (NO, NO2, NO3) implied in the found conservation things
nox_loc = (gc_spc_names .== "NO") .+ (gc_spc_names .== "NO2") .+ (gc_spc_names .== "NO3") 
cond(hcat(Γ_gc,nox_loc)'*hcat(Γ_gc,nox_loc))

# Application to POLLU
# from here: https://doi.org/10.1137/0915076
A_pollu = zeros(Int, 20, 25)
# Fill the matrix with coefficients
A_pollu[1, [1, 10, 14, 23, 24, 2, 3, 9, 11, 12, 22, 25]] .= [-1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1]
A_pollu[2, [2, 3, 9, 12, 1, 21]] .= [-1, -1, -1, -1, 1, 1]
A_pollu[3, [15, 1, 17, 19, 22]] .= [-1, 1, 1, 1, 1]
A_pollu[4, [2, 16, 17, 23, 15]] .= [-1, -1, -1, -1, 1]
A_pollu[5, [3, 4, 6, 7, 13, 20]] .= [-1, 2, 1, 1, 1, 1]
A_pollu[6, [6, 8, 14, 20, 3, 18]] .= [-1, -1, -1, -1, 1, 2]
A_pollu[7, [4, 5, 6, 13]] .= [-1, -1, -1, 1]
A_pollu[8, [4, 5, 6, 7]] .= [1, 1, 1, 1]
A_pollu[9, [7, 8]] .= [-1, -1]
A_pollu[10, [12, 7, 9]] .= [-1, 1, 1]
A_pollu[11, [9, 10, 8, 11]] .= [-1, -1, 1, 1]
A_pollu[12, [9]] .= [1]
A_pollu[13, [11, 10]] .= [-1, 1]
A_pollu[14, [13, 12]] .= [-1, 1]
A_pollu[15, [14]] .= [1]
A_pollu[16, [18, 19, 16]] .= [-1, -1, 1]
A_pollu[17, [20]] .= [-1]
A_pollu[18, [20]] .= [1]
A_pollu[19, [21, 22, 24, 23, 25]] .= [-1, -1, -1, 1, 1]
A_pollu[20, [25, 24]] .= [-1, 1]

pollu_spc_names = spc_names = ["[NO2]", "[NO]", "[O3P]", "[03]", "[HO2]", "[OH]", "[HCHO]", "[CO]", "[ALD]", "[MEO2]", "[C203]", "[CO2]", "[PAN]", "[CH30]", "[HNO3]", "[OID]", "[SO2]", "[SO4]", "[NO3]", "[N205]"]

Γ_pollu = conservationQR(A_pollu)
pollu_cycles = print_conserved_cycles(pollu_spc_names,Γ_pollu)