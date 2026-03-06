import math

def sigma(kappa1, alpha1, kappa2, alpha2):
    return (kappa1/math.sqrt(alpha1)) / (kappa2/math.sqrt(alpha2))

def PrandtlNum():
    return 0.7

def grashofNum(rho, g, mu, beta):
    numerator = g*(rho**2)*beta

    denominator = (mu**2)

    return numerator/denominator


## Solid
kappaS = 1600
cS = 500
rhoS = 8000

alphaS = kappaS/(cS*rhoS)

print(alphaS)
print("Alpha solid: {:f}".format(alphaS))

# Fluid
cF = 1000
muF = 7e-4
kappaF = cF*muF/PrandtlNum()
p = 1e5

beta = 0.01


rhoF = 1.0
alphaF = kappaF/(cF*rhoF)

print("rho fluid: {:f}".format(rhoF))
print("alpha fluid: {:f}".format(alphaF))


print("kappa fluid: {:f}".format(kappaF))
print("kS/kF: {:f}".format(kappaS/kappaF))
print("alphaS/alphaF: {:f}".format(alphaS/alphaF))


g =142857.142857/grashofNum(rhoF, 1, muF, beta)

print("g: {:f}".format(g))

print("Grashof number: {:f}".format(grashofNum(rhoF, g, muF, beta)))

print("Sigma: {:f}".format(sigma(kappaF, alphaF, kappaS, alphaS)))