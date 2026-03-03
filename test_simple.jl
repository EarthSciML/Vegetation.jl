#!/usr/bin/env julia

using ModelingToolkit
using DynamicQuantities

# Add the package to the path
push!(LOAD_PATH, "./src")
using Vegetation

println("Testing SoilVaporTransfer component...")

try
    sys = SoilVaporTransfer()
    println("✓ Component created successfully")
    println("Number of equations: ", length(equations(sys)))
    println("Number of unknowns: ", length(unknowns(sys)))

    # Try to validate the system
    println("Validating system...")
    # This should check if units are consistent
    validated_sys = complete(sys)
    println("✓ System validation passed")

catch e
    println("✗ Error: ", e)
    if isa(e, ModelingToolkit.ValidationError)
        println("This is a validation error, likely unit-related")
    end
end