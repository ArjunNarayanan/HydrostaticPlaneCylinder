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

function direct_core_strain_energy(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    A1c, A1s, A2s, B = HS.displacement_coefficients(
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

    A1c, A1s, A2s, B = HS.displacement_coefficients(
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
    A1c, A1s, A2s, B = HS.displacement_coefficients(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )

    return -V0c * (
        2 * (lambda + mu) * A1c^2 +
        2 * (lambda + mu) * A1c * B +
        (lambda - 2mu) / 2 * B^2 +
        2 * (lambda + mu) * A1c +
        lambda * B
    )
end

function expanded_core_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    K = bulk_modulus(lambda, mu)
    g = inner_radius / outer_radius
    C1 = (lambda - 2mu) / (lambda + 2mu)

    t1 = -(lambda / 3 - (lambda + mu) * C1 / 3) * theta0
    t2 = -(-lambda / 3 + (lambda + mu) * C1 / 3) * theta0 * g^2
    t3 =
        -(
            (lambda - 2mu) / 18 - (lambda + mu) * C1 / 9 +
            (lambda + mu) * C1^2 / 18
        ) * theta0^2
    t4 =
        -(
            -(lambda - 2mu) / 9 + 2 * (lambda + mu) * C1 / 9 -
            (lambda + mu) * C1^2 / 9
        ) *
        theta0^2 *
        g^2
    t5 =
        -(
            (lambda - 2mu) / 18 - (lambda + mu) / 9 * C1 +
            (lambda + mu) * C1^2 / 18
        ) *
        theta0^2 *
        g^4

    return V0c * (t1 + t2 + t3 + t4 + t5)
end

function simplified_core_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0c,
)
    C1 = (lambda - 2mu) / (lambda + 2mu)
    K = bulk_modulus(lambda, mu)
    C3 = mu * K / (lambda + 2mu)
    gamma = inner_radius / outer_radius

    return V0c*(-C3 * theta0 + C3 * theta0 * gamma^2 + C1 * C3 / 6 * theta0^2 -
    C1 * C3 / 3 * theta0^2 * gamma^2 + C1 * C3 / 6 * theta0^2 * gamma^4)
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
err = maximum_error(solverpotential, directpotential)

@test err < 10eps()

expandedpotential =
    expanded_core_potential.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0c,
    )

err = maximum_error(expandedpotential, solverpotential)
@test err < 10eps()

simplifiedpotential = simplified_core_potential.(inner_radius,outer_radius,lambda,mu,theta0,V0c)

err = maximum_error(simplifiedpotential,solverpotential)
@test err < 10eps()
