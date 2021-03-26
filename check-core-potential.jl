include("hydrostatic-solver.jl")
include("moduli-conversion.jl")
include("utilities.jl")
HS = HydrostaticSolver

function solver_core_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    solver = HS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    return HS.core_strain_energy(solver, V0c)
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

function direct_core_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    A1c, A1s, A2s, B = displacement_coefficients(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )

    cse =
        V0c * (
            2 * (lambda + mu) * A1c^2 +
            2 * lambda * A1c * B +
            (lambda + 2mu) * B^2 / 2
        )
    return cse
end

function solver_core_compression_work(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)

    solver = HS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    V = V0c * (1 + sum(HS.core_strain(solver)))
    srr = HS.core_stress(solver)[1]
    return V * srr
end

function direct_core_compression_work(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)

    A1c, A1s, A2s, B = displacement_coefficients(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )
    cw =
        V0c * (
            4 * (lambda + mu) * A1c^2 +
            4lambda * A1c * B +
            2mu * A1c * B +
            lambda * B^2 +
            2(lambda + mu) * A1c +
            lambda * B
        )
    return cw
end

function solver_core_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    solver = HS.CylindricalSolver(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        lambda,
        mu,
        theta0,
    )
    cse = HS.core_strain_energy(solver, V0c)
    ccw = HS.core_compression_work(solver, V0c)
    return cse - ccw
end

function direct_core_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    A1c, A1s, A2s, B = displacement_coefficients(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )

    return V0c * (
        -2 * (lambda + mu) * A1c^2 - 2 * (lambda + mu) * A1c * B -
        (lambda - 2mu) / 2 * B^2 - 2 * (lambda + mu) * A1c - lambda * B
    )
end

function maximum_error(u, v)
    return maximum(abs.(u - v))
end

K = 247.0
mu = 126.0
lambda = lame_lambda(K, mu)
theta0 = 0.067
rhoc = 3.68e3
rhos = 3.93e3

V0c = inv(rhoc)
V0s = inv(rhos)

outer_radius = 1.0
dx = outer_radius / 1e3
inner_radius = dx:dx:outer_radius

solvercse =
    solver_core_strain_energy.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )
directcse =
    direct_core_strain_energy.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )

err = maximum_error(solvercse, directcse)

using Test
@test err < 10eps()

solvercw =
    solver_core_compression_work.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )
directcw =
    direct_core_compression_work.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )

err = maximum_error(solvercw, directcw)
@test err < 10eps()

solverpotential =
    solver_core_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0c)
directpotential =
    direct_core_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0c)

err = maximum_error(solverpotential,directpotential)

@test err < 10eps()
