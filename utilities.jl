function pressure(stress)
    return -1.0 / 3.0 * sum(stress)
end

function deviatoric_stress(stress)
    p = pressure(stress)
    return stress .+ p
end

function deviatoric_strain(strain)
    return strain .- 1.0/3.0*sum(strain)
end

function second_invariant(strain)
    return 0.5*sum(strain.^2)
end
