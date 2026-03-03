#!/usr/bin/env julia

using ModelingToolkit
using DynamicQuantities
using Vegetation

println("Testing basic component creation...")
try
    sys = SoilVaporTransfer()
    println("✓ SoilVaporTransfer component created successfully")
    println("Number of equations: ", length(equations(sys)))
    println("Number of unknowns: ", length(unknowns(sys)))
catch e
    println("✗ Error creating SoilVaporTransfer component:")
    println(e)
end

println("\nTesting PDE system creation (no vapor)...")
try
    pde_prel = SoilVaporTransferPDE(0.5, 3600.0; include_vapor = false)
    println("✓ PDE system (no vapor) created successfully")
catch e
    println("✗ Error creating PDE system (no vapor):")
    println(e)
end

println("\nTesting PDE system creation (with vapor)...")
try
    pde_simp = SoilVaporTransferPDE(0.5, 3600.0; include_vapor = true)
    println("✓ PDE system (with vapor) created successfully")
catch e
    println("✗ Error creating PDE system (with vapor):")
    println(e)
end

println("\nDone.")