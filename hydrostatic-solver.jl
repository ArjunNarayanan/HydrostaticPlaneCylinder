module HydrostaticSolver

include("moduli-conversion.jl")
include("utilities.jl")

function coefficient_matrix(inner_radius, outer_radius, ls, ms, lc, mc)
    R = inner_radius
    b = outer_radius
    R2 = R^2
    b2 = b^2

    row1 = [R, -R, -1 / R, 0]
    row2 = [2(lc + mc), -2(ls + ms), 2ms / R2, lc - ls]
    row3 = [0, 2(ls + ms), -2ms / b2, ls]
    row4 =
        [2lc * R2, 2ls * (b2 - R2), 0, (lc + 2mc) * R2 + (ls + 2ms) * (b2 - R2)]

    op = vcat(row1', row2', row3', row4')
end

function coefficient_rhs(inner_radius, outer_radius, ls, ms, theta0)
    K = bulk_modulus(ls, ms)
    v = [
        0,
        -K * theta0,
        K * theta0,
        K * theta0 * (outer_radius^2 - inner_radius^2),
    ]
    return v
end

struct CylindricalSolver
    inner_radius::Any
    outer_radius::Any
    ls::Any
    ms::Any
    lc::Any
    mc::Any
    theta0::Any
    A1c::Any
    A1s::Any
    A2s::Any
    B::Any
    function CylindricalSolver(
        inner_radius,
        outer_radius,
        ls,
        ms,
        lc,
        mc,
        theta0,
    )
        m = coefficient_matrix(inner_radius, outer_radius, ls, ms, lc, mc)
        r = coefficient_rhs(inner_radius, outer_radius, ls, ms, theta0)
        A1c, A1s, A2s, B = m \ r
        new(
            inner_radius,
            outer_radius,
            ls,
            ms,
            lc,
            mc,
            theta0,
            A1c,
            A1s,
            A2s,
            B,
        )
    end
end

function core_strain(solver::CylindricalSolver)
    return [solver.A1c, solver.A1c, solver.B]
end

function shell_strain(solver::CylindricalSolver, r)
    A1 = solver.A1s
    A2 = solver.A2s
    t0 = solver.theta0
    B = solver.B

    err = A1 - A2 / r^2 - t0 / 3
    ett = A1 + A2 / r^2 - t0 / 3
    ezz = B - t0 / 3
    return [err, ett, ezz]
end

function core_stress(solver::CylindricalSolver)
    l = solver.lc
    m = solver.mc
    A1 = solver.A1c
    B = solver.B

    srr = 2 * (l + m) * A1 + l * B
    stt = srr
    szz = 2 * l * A1 + (l + 2 * m) * B

    return [srr, stt, szz]
end

function shell_stress(solver::CylindricalSolver, r)
    l = solver.ls
    m = solver.ms
    t0 = solver.theta0
    A1 = solver.A1s
    A2 = solver.A2s
    B = solver.B

    KT = bulk_modulus(l, m) * t0

    srr = 2 * (l + m) * A1 - 2 * m * A2 / r^2 + l * B - KT
    stt = 2 * (l + m) * A1 + 2 * m * A2 / r^2 + l * B - KT
    szz = 2 * l * A1 + (l + 2 * m) * B - KT

    return [srr, stt, szz]
end

function core_strain_energy(solver::CylindricalSolver, V0c)
    stress = core_stress(solver)
    strain = core_strain(solver)
    return 0.5 * V0c * sum(stress .* strain)
end

function shell_strain_energy(solver::CylindricalSolver, r, V0s)
    stress = shell_stress(solver, r)
    strain = shell_strain(solver, r)

    return 0.5 * V0s * sum(stress .* strain)
end

function core_compression_work(solver::CylindricalSolver, V0)
    V = V0 * (1.0 + sum(core_strain(solver)))
    srr = core_stress(solver)[1]
    return V * srr
end

function shell_compression_work(solver::CylindricalSolver, r, V0)
    V = V0 * (1.0 + sum(shell_strain(solver, r)))
    srr = shell_stress(solver, r)[1]
    return V * srr
end

function displacement_coefficients(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
)
    K = bulk_modulus(lambda, mu)
    C1 = (lambda - 2mu) / (lambda + 2mu)
    C2 = K / (2 * (lambda + 2mu))
    g = inner_radius / outer_radius

    A1c = C1 * theta0 / 6 * (g^2 - 1)
    A1s = C1 * theta0 / 6 * (g^2 + 2 / C1)
    A2s = -C2 * theta0 * inner_radius^2
    B = theta0 / 3 * (1 - g^2)
    return A1c, A1s, A2s, B
end

end
