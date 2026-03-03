module Vegetation

using EarthSciMLBase
using ModelingToolkit
using ModelingToolkit: t, D
using DynamicQuantities
using DocStringExtensions

include("landis_biomass.jl")
export LANDISBiomass

include("stage_prognosis.jl")
export StagePrognosis, StagePrognosisHCB

include("maize_root_growth.jl")
export MaizeRootGrowth

end
