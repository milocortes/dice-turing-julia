using Turing
using Statistics: mean, std
using Optim
#using Random:seed!
using StatsPlots
#seed!(123)

#set
t = Vector(1:100)
NT = length(t)

# PARÁMETROS

#****** Disponibilidad de combustibles fósiles
fosslim = 6000 # Extracción acumulativa máxima de combustibles fósiles (GtC)  /6000/

#****** Paso de tiempo
Δ = 5 # Años por periodo /5/

#****** If optimal control
#        ifopt  = 0 #  Indicator where optimized is 1 and base is 0      /0/

#****** Preferencias
α = 1.45 # Elasticidad de la utilidad marginal al consumo /1.45/
ρ = 0.015 # Tasa de descuento por año /0.015/

#****** Población y tecnología
γ = 0.300   # Elasticidad de sustitución del capital en la función de producción /0.300/
pop0 = 7403 # Población inicial mundial de 2015 (millones) /7403/
𝓁_g = 0.134 # Tasa de crecimiento para calibrar la proyección de la población de 2050 /0.134/
Lₐ = 11500 # Población asintótica (millones) /11500/
δₖ = 0.100 # Tasa de depreciación del capital (por año) /0.100/
q0 = 105.5 # Producto bruto inicial mundial 2015 (trill 2010 USD)/105.5/
k0 = 223   # Valor del capital inicial de 2015 (trill 2010 USD)/223/
a0 = 5.115 # Nivel inicial de la productividad total de los factores /5.115/
gₐ = 0.076 # Tasa de crecimiento inicial para la PTF por 5 años /0.076/
δₐ = 0.005 # Tasa de caída de la PTF por 5 años /0.005/

#****** Parámetros de emisiones
gσ = -0.0152 # Crecimiento inicial de sigma (por año) /-0.0152/
δσ = -0.001 # Tasa de caída de descarbonización (por periodo) /-0.001/
θ₂ = 2.6 # Emisiones de carbono de la tierra en 2015 (GtC02 por año) /2.6/
δEL = 0.115 # Tasa de caída de las emisiones de la tierra (por periodo) /0.115/
e0 = 35.85 # Emisiones industriales de 2015 (GtC02 por año) /35.85/
miu0 = 0.03 # Tasa inicial de control de las emisiones para el caso base 2015 /0.03/

#****** Ciclo del carbón

#-- Condiciones iniciales
mat0 = 851 #  Initial Concentration in atmosphere 2015 (GtC)       /851  /
mu0  = 460 #  Initial Concentration in upper strata 2015 (GtC)     /460  /
ml0  = 1740 #  Initial Concentration in lower strata 2015 (GtC)    /1740 /
mateq = 588 # mateq Equilibrium concentration atmosphere  (GtC)    /588  /
mueq  = 360 # mueq Equilibrium concentration in upper strata (GtC) /360  /
mleq = 1720 # mleq Equilibrium concentration in lower strata (GtC) /1720 /

#-- Parámetros de flujo
𝜁₁₂  = .12 #    Carbon cycle transition matrix                     /.12  /
𝜁₂₃  = 0.007 #   Carbon cycle transition matrix                    /0.007/

#-- These are for declaration and are defined later
𝜁₁₁  = 1 - 𝜁₁₂   # Carbon cycle transition matrix
𝜁₂₁  = 𝜁₁₂*mateq / mueq   # Carbon cycle transition matrix
𝜁₂₂  = 1 - 𝜁₂₁ - 𝜁₂₃ # Carbon cycle transition matrix
𝜁₃₂  = 𝜁₂₃*mueq / mleq # Carbon cycle transition matrix
𝜁₃₃  = 1 - 𝜁₃₂  # Carbon cycle transition matrix
σ₀  = e0/(q0*(1-miu0))  # Carbon intensity 2010 (kgCO2 per output 2005 USD 2010)

#****** Parámetros del modelo climático

t2xco2  = 3.1 # Equilibrium temp impact (oC per doubling CO2)    / 3.1 /
fex0  = 0.5 #   2015 forcings of non-CO2 GHG (Wm-2)              / 0.5 /
fex1  = 1.0 #   2100 forcings of non-CO2 GHG (Wm-2)              / 1.0 /
tocean0  = .0068 # Initial lower stratum temp change (C from 1900) /.0068/
tatm0  = 0.85 #  Initial atmospheric temp change (C from 1900)    /0.85/
c1  = 0.1005 #     Climate equation coefficient for upper level  /0.1005/
c3  = 0.088 #     Transfer coefficient upper to lower stratum    /0.088/
c4  = 0.025 #     Transfer coefficient for lower level           /0.025/
F₂ₓ  = 3.6813 # Forcings of equilibrium CO2 doubling (Wm-2)   /3.6813 /

#****** Parámetros de daño climático
ψ₁₀ = 0 #     Initial damage intercept                         /0   /
ψ₁  = 0 #      Damage intercept                                 /0   /
ψ₂  = 0.00236 #      Damage quadratic term                     /0.00236/
ψ₂₀ = ψ₂ #     Initial damage quadratic term
ψ₃  = 2.00 #      Damage exponent                              /2.00   /

#** Abatement cost
expcost2 = 2.6 # Exponent of control cost function             / 2.6  /
pback  = 550 #   Cost of backstop 2010$ per tCO2 2015          / 550  /
gback  = .025 #   Initial cost decline backstop cost per period / .025/
limmiu  = 1.2 #  Upper limit on control rate after 2150        / 1.2 /
tnopol  = 45 #  Period before which no emissions controls base  / 45   /
cprice0  = 2 # Initial base carbon price (2010$ per tCO2)      / 2    /
gcprice  = .02 # Growth rate of base carbon price per year     /.02  /
#** Scaling and inessential parameters
#* Note that these are unnecessary for the calculations
#* They ensure that MU of first period's consumption =1 and PV cons = PV utilty
scale1  = 0.0302455265681763 #    Multiplicative scaling coefficient           /0.0302455265681763 /
scale2  = -10993.704 #    Additive scaling coefficient       /-10993.704/;
#* Program control variables
lam = F₂ₓ/ t2xco2


# l(t) Nivel de la población y del trabajo
function initPoblacion(pop0,Lₐ,𝓁_g,NT)
    l = zeros(NT)
    l[1] = pop0
    for i in 2:NT
        l[i] = l[i-1] * (Lₐ/l[i-1])^𝓁_g
    end

    return l
end

# al(t) Nivel de la productividad total de los factores
function initPTF(gₐ,δₐ,a0,t)
    ga = gₐ*ℯ.^(-δₐ .*5 .*(t.-1))
    NT = length(t)
    al = zeros(NT)
    al[1] = a0

    for i in 2:NT
        al[i] = al[i-1]/(1-ga[i-1])
    end

    return al,ga
end

function initGσ(gσ,NT,δσ,Δ)
    gsig = zeros(NT)
    gsig[1] = gσ

    for i in 2:NT
        gsig[i] = gsig[i-1]*((1+δσ)^Δ)
    end

    return gsig
end

function initσ(σ₀,NT,gsig,Δ)
    σ = zeros(NT)
    σ[1] = σ₀

    for i in 2:NT
        σ[i] = σ[i-1]*ℯ^(gsig[i-1]*Δ)
    end

    return σ
end

function initCarbonTree(θ₂,δEL,t)
    etree = θ₂*(1-δEL).^(t.-1)
    NT = length(t)
    cumetree =  zeros(NT)
    cumetree[1] = 100
    for i in 2:NT
        cumetree[i] = cumetree[i-1] + etree[i-1]*(5/3.666)
    end

    return etree,cumetree
end

l = initPoblacion(pop0,Lₐ,𝓁_g,NT)
al,ga = initPTF(gₐ,δₐ,a0,t)
gsig = initGσ(gσ,NT,δσ,Δ)
sigma = initσ(σ₀,NT,gsig,Δ)
etree,cumetree = initCarbonTree(θ₂,δEL,t)
pbacktime = pback * (1-gback) .^(t.-1)
cost1 = pbacktime .* sigma / expcost2 / 1000
rr =  1 ./((1+ρ).^(Δ .*(t .-1 )))
forcoth = fill(fex0, NT)
forcoth[1:17] = forcoth[1:17] .+ (1/17) .*(fex1-fex0)*(t[1:17] .-1)
forcoth[18:NT] = forcoth[18:NT] .+ (fex1-fex0)

optlrsav = (δₖ + .004)/(δₖ + .004*α + ρ)*γ
cpricebase = cprice0 .*(1 +gcprice) .^(5 .*(t .-1))


## Definimos arreglos
K = zeros(NT)
K[1] = k0
CCA = zeros(NT)
CCA[1] = 400
MAT = zeros(NT)
MAT[1] = mat0
ML = zeros(NT)
ML[1] <- ml0
MU = zeros(NT)
MU[1] = mu0
TATM = zeros(NT)
TATM[1] = tatm0
TOCEAN = zeros(NT)
TOCEAN[1] = tocean0

YGROSS = zeros(NT)
EIND = zeros(NT)
CCATOT = zeros(NT)
FORC = zeros(NT)
DAMFRAC = zeros(NT)
DAMAGES = zeros(NT)
ABATECOST = zeros(NT)
MCABATE = zeros(NT)
YNET = zeros(NT)
Y = zeros(NT)
I = zeros(NT)
C = zeros(NT)
CPC = zeros(NT)
PERIODU = zeros(NT)
CEMUTOTPER = zeros(NT)
E = zeros(NT)

# * Control rate limits
MIU = zeros(NT)
dMIU = zeros(NT)
miu_lo = fill(0.1, NT)
miu_up = fill(0.999999, NT)

S = zeros(NT)
s_lo = fill(1e-1, NT)
s_up = fill(0.9, NT)
lag10 =  t .> NT -10
s_up[lag10] .= optlrsav
tatm_sampling = zeros(NT)

# Definimos el modelo de Turing
@model max_DICE(NT) =
begin

    # Priors
    S_tau = 0.05
    MIU_tau = 0.05
    dMIU[1] ~ Uniform(miu_lo[1], miu_up[1])
    MIU[1] = dMIU[1]

    S[1] ~ Uniform(s_lo[1], s_up[1])

    # Random walk priors
    for i in 2:NT
        #dMIU[i] ~ Normal(0, MIU_tau/rr[i])
        #MIU[i] = min(max(MIU[i-1] + dMIU[i], miu_lo[i]), miu_up[i])
        dMIU[i] ~ Uniform(0 , (MIU_tau/rr[i]))
        MIU[i] = min(max(MIU[i-1]+ dMIU[i],miu_lo[i]),miu_up[i])
        S[i] ~ Truncated(Normal(S[i-1], (rr[i]/S_tau)),s_lo[i], s_up[i])
    end

    # Likelihood
    for i in 1:NT-1
        K[i+1] = max(0,((1-δₖ)^Δ)*K[i] + Δ*I[i])
        YGROSS[i] =  al[i] * ((l[i]/1000)^(1-γ)) * K[i]^γ
        EIND[i] = sigma[i] * YGROSS[i] * (1 - MIU[i])
        E[i] = EIND[i] + etree[i] #+ dE
        CCA[i+1] = CCA[i] + EIND[i] * 5 / 3.666
        CCATOT[i] = CCA[i] + cumetree[i]
        MAT[i+1] = MAT[i]*𝜁₁₁ + MU[i]*𝜁₂₁ + E[i] * 5 / 3.666
        ML[i+1] = ML[i]*𝜁₃₃  + MU[i]*𝜁₂₃
        MU[i+1] = MAT[i]*𝜁₁₂ + MU[i]*𝜁₂₂ + ML[i]*𝜁₃₂
        FORC[i] = F₂ₓ * log(MAT[i]/588.000)/log(2) + forcoth[i]
        TATM[i+1] = TATM[i] + c1 * (FORC[i+1] - (F₂ₓ/t2xco2) * TATM[i] - c3 * (TATM[i] - TOCEAN[i]))
        TOCEAN[i+1] = TOCEAN[i] + c4*(TATM[i]-TOCEAN[i])
        DAMFRAC[i] = ψ₁*TATM[i] + ψ₂*TATM[i]^ψ₃
        DAMAGES[i] = YGROSS[i] * DAMFRAC[i]
        ABATECOST[i] = YGROSS[i] * cost1[i] * MIU[i]^expcost2
        MCABATE[i] = pbacktime[i] * MIU[i]^(expcost2-1)
        YNET[i] = YGROSS[i] * (1 - DAMFRAC[i])
        Y[i] = YNET[i] - ABATECOST[i]
        I[i] = S[i] * Y[i]
        C[i] = Y[i] - I[i]
        CPC[i] = 1000 * C[i] / l[i]
        PERIODU[i] = ((C[i]*1000/l[i])^(1-α) - 1) / (1 - α) - 1
        CEMUTOTPER[i] = PERIODU[i] * l[i] * rr[i]
    end;

    # Calculamos para NT
    YGROSS[NT] =  al[NT] * ((l[NT]/1000)^(1-γ)) * K[NT]^γ
    EIND[NT] = sigma[NT] * YGROSS[NT] * (1 - MIU[NT])
    E[NT] = EIND[NT] + etree[NT] #+ dE
    CCATOT[NT] = CCA[NT] + cumetree[NT]
    FORC[NT] = F₂ₓ * log(MAT[NT]/588.000)/log(2) + forcoth[NT]
    DAMFRAC[NT] = ψ₁*TATM[NT] + ψ₂*TATM[NT]^ψ₃
    DAMAGES[NT] = YGROSS[NT] * DAMFRAC[NT]
    ABATECOST[NT] = YGROSS[NT] * cost1[NT] * MIU[NT]^expcost2
    MCABATE[NT] = pbacktime[NT] * MIU[NT]^(expcost2-1)
    YNET[NT] = YGROSS[NT] * (1 - DAMFRAC[NT])
    Y[NT] = YNET[NT] - ABATECOST[NT]
    I[NT] = S[NT] * Y[NT]
    C[NT] = Y[NT] - I[NT]
    CPC[NT] = 1000 * C[NT] / l[NT]
    PERIODU[NT] = ((C[NT]*1000/l[NT])^(1-α) - 1) / (1 - α) - 1
    CEMUTOTPER[NT] = PERIODU[NT] * l[NT] * rr[NT]



    # Calculamos la utilidad del periodo
    UTILITY = Δ * scale1 * sum(CEMUTOTPER) + scale2
    OBJ ~ Normal(UTILITY, 100)
    println(OBJ)

    return (;TATM)

end;

# I use optim for the maximum Likelihood estimation
model_optim = max_DICE(NT)
obj_max = 0
mle_estimate = 0
while obj_max < 4000
    #mle_estimate = optimize(model_optim,MAP(), NelderMead())
    mle_estimate = optimize(model_optim,MLE(), NelderMead(), Optim.Options(iterations=10_000, allow_f_increases=true))
    obj_max = mle_estimate.values[:OBJ]
end

## TATM
plot([2000+ i*5 for i in 1:100],TATM, title = "Increase temperature of the atmosphere (TATM)")
xlabel!("Year")
ylabel!("Degrees C from 1900")
savefig("tatm.png")

## DAMFRAC
plot([2000+ i*5 for i in 1:100],DAMFRAC, title = "Damages as fraction of gross output")
xlabel!("Year")
ylabel!("% of gross output")
savefig("damfrac.png")

## I sampling the chain with the optimization vector parameters
mc_optim = sample(model_optim, SMC(), MCMCThreads(), 10000, 4,init_theta = mle_estimate.values.array)

## TATM
plot([2000+ i*5 for i in 1:100],TATM, title = "Increase temperature of the atmosphere (TATM) (MCMC)")
xlabel!("Year")
ylabel!("Degrees C from 1900")
savefig("tatm_mcmc.png")

## DAMFRAC
plot([2000+ i*5 for i in 1:100],DAMFRAC, title = "Damages as fraction of gross output (MCMC)")
xlabel!("Year")
ylabel!("% of gross output")
savefig("damfrac_mcmc.png")
