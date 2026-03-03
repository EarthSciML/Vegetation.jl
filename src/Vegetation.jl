module Vegetation

using EarthSciMLBase
using ModelingToolkit
using ModelingToolkit: t, D, @register_symbolic
# Symbolics is needed for @register_symbolic macro expansion.
# Access it from loaded modules (transitive dep of ModelingToolkit) to avoid
# direct dependency version conflicts.
const Symbolics = Base.loaded_modules[
    Base.PkgId(
        Base.UUID("0c5d862f-8b57-4792-8d23-62f2024744c7"), "Symbolics"
    ),
]
using DynamicQuantities
using DocStringExtensions
using DomainSets: Interval
import RuntimeGeneratedFunctions

RuntimeGeneratedFunctions.init(@__MODULE__) # Needed even though we don't use it directly.

include("landis_biomass.jl")
export LANDISBiomass

include("stage_prognosis.jl")
export StagePrognosis, StagePrognosisHCB

include("vapor_transfer.jl")

end
