from sympy.physics.units import *
from sympy import *

# Rounding:
import decimal
from decimal import Decimal as DX
def iso_round(obj, pv, rounding=decimal.ROUND_HALF_EVEN):
    import sympy
    """
    Rounding acc. to DIN EN ISO 80000-1:2013-08
    place value = Rundestellenwert
    """
    assert pv in set([
        # place value   #  round to:
        100,            #  3rd last digit before decimal
        10,             #  2nd last
        1,              #  last
        0.1,            #  1st digit after decimal
        0.01,           #  2nd
        0.001,          #  3rd
        0.0001,         #  4th
        0.00001,        #  5th
        0.000001,       #  6th
        0.0000001,      #  7th
        0.00000001,     #  8th
        0.000000001,    #  9th
        0.0000000001,   # 10th
        ])
    try:
        tmp = DX(str(float(obj)))
        obj = tmp.quantize(DX(str(pv)), rounding=rounding)
    except:
        for i in range(len(obj)):
            tmp = DX(str(float(obj[i])))
            obj[i] = tmp.quantize(DX(str(pv)), rounding=rounding)
    return obj

# LateX:
kwargs = {}
kwargs["mat_str"] = "bmatrix"
kwargs["mat_delim"] = ""
# kwargs["symbol_names"] = {FB: "F^{\mathsf B}", }

# Units:
(k, M, G ) = ( 10**3, 10**6, 10**9 )
(mm, cm, deg) = ( m/1000, m/100, pi/180)
Newton = kg*m/s**2
Pa     = Newton/m**2
MPa    = M*Pa
GPa    = G*Pa
kN     = k*Newton

# ---

c, mass = var("c, mass", positive=True)

K = Matrix([
[2*c,  -c],
[-c,  c],
])

M = Matrix([
[3*mass/2, 0],
[0, mass],
])


EI = 210 *GPa * 5100 *cm**4
l = S(45)/10 *m
sub_list = [
    (mass, 1000*kg    ),
    (c, 24*EI/l**3 ),
    ]

# ξ = λ²
xi = var("xi")
w = var("omega")

A = K + xi*M

pprint("\nCharacteristic equation:")
eq = Eq(det(A))
pprint(eq)

sol_xi = solve(eq,xi)

w2, w3 = var("w2, w3")
w = Matrix([w2, w3])

zero = Matrix([0,0])
for i in range(len(sol_xi)):
    pprint("\n\nEigenvalue:")
    xii = sol_xi[i]
    pprint(xii)
    Ai = A.subs(xi,xii)
    eq = Eq(Ai*w,zero)
    sol = solve(eq, w2)
    pprint("\nEigenvector:")
    pprint(sol)

# Omega 1:
pprint("\nw1 / (1/s):")
w1 = sqrt(2*c/mass)
tmp = w1.subs(sub_list)
tmp /= (1/s)
tmp = iso_round(tmp, 0.1)
pprint(tmp)


w = w1.subs(sub_list)
w_in_Hz = w / (1/s)

pprint("\nPeriod T1 / s:")
T = 2*pi/w
T_in_s = T / s
T_in_s = N(T_in_s,20)
T_in_s = float(T_in_s)
tmp = T_in_s
tmp = iso_round(tmp,0.001)
pprint(tmp)

from pylab import *
from numpy import linspace
t_in_s = linspace(0,T_in_s,100)
wt =  w_in_Hz * t_in_s
cos_wt = array([cos(float(x)) for x in wt])
w2_in_mm = -10 * cos_wt
w3_in_mm = +10 * cos_wt


plt.axis()
plt.grid()
plt.plot(t_in_s, w2_in_mm, "b-", label=r"$w2\,\, /  \,\, \mathrm{mm}$")
plt.plot(t_in_s, w3_in_mm, "r--", label=r"$w3\,\, /  \,\, \mathrm{mm}$")
plt.xlabel(r"$t  \,\, /  \,\, \mathrm{s}$")
plt.legend()
plt.savefig('2dofs_motion.svg', transparent=True)
plt.show()

# Characteristic equation:
#    2                ⎛      3⋅mass⋅ξ⎞
# - c  + (c + mass⋅ξ)⋅⎜2⋅c + ────────⎟ = 0
#                     ⎝         2    ⎠
#
#
# Eigenvalue:
# -2⋅c
# ─────
#  mass
#
# Eigenvector:
# {w₂: -w₃}
#
#
# Eigenvalue:
#  -c
# ──────
# 3⋅mass
#
# Eigenvector:
# ⎧    2⋅w₃⎫
# ⎨w₂: ────⎬
# ⎩     3  ⎭
#
# w1 / (1/s):
# 75.1
#
# Period T1 / s:
# 0.084
