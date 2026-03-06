import math

def sigma(kappa1, alpha1, kappa2, alpha2):
    return (kappa1/math.sqrt(alpha1)) / (kappa2/math.sqrt(alpha2))

def PrandtlNum():
    return 7.0

def grashofNum(g, nu, beta):
    numerator = g*beta

    denominator = (nu**2)

    return numerator/denominator


def timeRequired(oldAlpha, newAlpha):
    return (oldAlpha/newAlpha)*0.07


## Solid
kappaS = 1600
cS = 50
rhoS = 750

alphaS = kappaS/(cS*rhoS)

print(alphaS)
print("Alpha solid: {:f}".format(alphaS))

# Fluid
cF = 175
muF = 0.8
kappaF = cF*muF/PrandtlNum()
p = 1e5


rhoF = 240

nuF = muF/rhoF
beta = 0.01



M = (8.314*273.15/p)

alphaF = kappaF/(cF*rhoF)

print("Molar mass {:f}".format(M))
print("rho fluid: {:f}".format(rhoF))
print("alpha fluid: {:f}".format(alphaF))

print("Dynamic Viscosity: {:f}".format(muF))
print("Kinematic Viscosity: {:f}".format(nuF))

print("kappa fluid: {:f}".format(kappaF))
print("kS/kF: {:f}".format(kappaS/kappaF))
print("alphaS/alphaF: {:f}".format(alphaS/alphaF))


g =10000/grashofNum(1, nuF, beta)

print("g: {:f}".format(g))

print("Grashof number: {:f}".format(grashofNum(g, nuF, beta)))

print("Sigma: {:f}".format(sigma(kappaF, alphaF, kappaS, alphaS)))