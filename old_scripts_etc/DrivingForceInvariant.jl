# does an invariant arise from linear dependency in the driving force?
using LinearAlgebra
using RowEchelon
include("/Users/psturm/Desktop/Julia Photochemical Model/assign2_rk.jl")
tempk = 298
press = 1
rk = zeros(Float64,10)
rk = getrk!(tempk,press,rk)

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

# r = DPx
D = Diagonal(rk)
P = zeros(Float64,10,9)
[P[i,i] = 1 for i=1:4]
[P[i+1,i] = 1 for i=4:9]

BDP = B*D*P
invariants = nullspace(BDP')

HC = [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0]
HN = [0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0]
CQ3_snapped = [-6, 5, -1, -3, -9, -6, -3, -6, -4, 3, 2.213115]
cq_charlie = [-2,3,1,-2,-3,-2,-1,-2,0,0, 4*rk[4]/rk[5] - 2]


all_invariants = hcat(invariants,HC,HN)
# all_invariants = hcat(invariants,CQ3_snapped)
cond(all_invariants'*all_invariants)

# plot(X[2:10000,5:end]*invariants[:,3])