#INCLUDE atoms.kpp

{ ========================= }
{ UCI Chem species          }
{ ========================= }
#DEFVAR
O3        = 3O;                     { Ozone }
OH        = O + H;
HO2       = 2O + H;
H2O2      = 2O + 2H;
CH2O      = C + 2H + O;
CH3O2     = C + 3H + 2O;
CH3OOH    = C + 4H + 2O;
NO        = N + O;
NO2       = N + 2O;
NO3       = N + 3O;
N2O5      = 2N + 5O;
HNO3      = N + H + 3O;
HO2NO2    = H + N + 4O;
PAN       = 2C + 3H + N + 5O;       { CH3CO3NO2 }
CO        = C + O;
C2H6      = 2C + 6H;
C3H8      = 3C + 8H;
C2H4      = 2C + 4H;
ROHO2     = 2C + 5H + 3O;           { C2H5O3 }
CH3COCH3  = 3C + 6H + O;
C2H5O2    = 2C + 5H + 2O;
C2H5OOH   = 2C + 6H + 2O;
CH3CHO    = 2C + 4H + O;
CH3CO3    = 2C + 3H + 3O;
ISOP      = 5C + 8H;
ISOP_VBS  = 5C + 8H;
ISOPO2    = 6C + 9H + 3O;           { HOCH2COOCH3CHCH2 }
C10H16    = 10C + 16H;
MVKMACR  = 8C + 12H + 2O;           { CH2CHCOCH3 + CH2CCH3CHO }
MVKO2    = 4C + 7H + 4O;
E90      = 3O;
N2OLNZ   = 2N;
NOYLNZ   = N;
CH4LNZ   = C + 4H;
H2OLNZ   = 2H + O;
DMS       = 2C + 6H + S;
SO2       = S + 2O;
H2SO4     = 2H + S + 4O;
SOAG0    = 15C + 38H + 2O;
SOAG15   = 15C + 38H + 2O;
SOAG24   = 15C + 38H + 2O;
SOAG31   = 15C + 38H + 2O;
SOAG32   = 15C + 38H + 2O;
SOAG33   = 15C + 38H + 2O;
SOAG34   = 15C + 38H + 2O;
SOAG35   = 15C + 38H + 2O;
so4_a1   = N + 4H + S + 4O;          { NH4HSO4 }
so4_a2   = N + 4H + S + 4O;
so4_a3   = N + 4H + S + 4O;
so4_a5   = N + 4H + S + 4O;
pom_a1   = C;
{ pom_a2 does not exist }
pom_a3   = C;
pom_a4   = C;
soa_a1   = 15C + 38H + 2O;
soa_a2   = 15C + 38H + 2O;
soa_a3   = 15C + 38H + 2O;
bc_a1    = C;
{ bc_a2 does not exist }
bc_a3    = C;
bc_a4    = C;
dst_a1   = Al + Si + 5O;
{ dst_a2 does not exist }
dst_a3   = Al + Si + 5O;
{ dst_a4 does not exist }
ncl_a1   = Na + Cl;
ncl_a2   = Na + Cl;
ncl_a3   = Na + Cl;
mom_a1   = 8520C + 11360H + 8520O;
mom_a2   = 8520C + 11360H + 8520O;
mom_a3   = 8520C + 11360H + 8520O;
mom_a4   = 8520C + 11360H + 8520O;
num_a1   = H;
num_a2   = H;
num_a3   = H;
num_a4   = H;
num_a5   = H;




#DEFFIX
M          = IGNORE;               { Atmospheric generic molecule }
N2         = 2N;
O2         = 2O;
H2O        = 2H + O;
H2         = 2H;
CH4        = C + 4H;
prsd_O3   = 3O;
prsd_NO3  = N + 3O;
prsd_OH   = O + H;
H         = H; {added by Obin Jan 2026}
CH3CO     = 2C + 3H + O;        {added by Obin Jan 2026}
