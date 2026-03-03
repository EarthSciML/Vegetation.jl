module Vegetation

using EarthSciMLBase
using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using DocStringExtensions
using DomainSets: Interval

include("landis_biomass.jl")
export LANDISBiomass

include("stage_prognosis.jl")
export StagePrognosis, StagePrognosisHCB

include("residue_mulch.jl")
export ResidueMulchDecomposition, MulchRadiationAttenuation, MulchWindProfile,
    MulchHeatVaporFluxes, MulchWaterCharacteristic,
    MulchHeatWaterTransfer, MulchHeatWaterPDE, MulchSurfaceRunoffPDE

end
