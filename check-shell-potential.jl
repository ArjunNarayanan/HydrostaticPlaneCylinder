using Test
include("hydrostatic-solver.jl")
include("moduli-conversion.jl")
include("utilities.jl")
HS = HydrostaticSolver

function solver_shell_strain_energy(
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
    sse = HS.shell_strain_energy(solver, inner_radius)
    return sse
end

function direct_shell_strain_energy(
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
    K = bulk_modulus(lambda, mu)

    sse = (
        2 * (lambda + mu) * (A1s)^2 +
        2mu * (A2s / inner_radius^2)^2 +
        (lambda + 2mu) * B^2 / 2 +
        2lambda * A1s * B - K * theta0 * B - 2K * theta0 * A1s +
        K * theta0^2 / 2
    )
    return sse
end

function solver_shell_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0s,
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
    sse = HS.shell_strain_energy(solver, inner_radius, V0s)
    scw = HS.shell_compression_work(solver, inner_radius, V0s)
    return sse - scw
end

function direct_shell_potential(
    inner_radius,
    outer_radius,
    lambda,
    mu,
    theta0,
    V0s,
)
    A1c, A1s, A2s, B = HS.displacement_coefficients(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
    )

    pd =
        -2 * (lambda + mu) * A1s^2 + 2mu * (A2s / inner_radius^2)^2 -
        (lambda - 2mu) / 2 * B^2 + 4mu * A1s * A2s / inner_radius^2 -
        2 * (lambda + mu) * A1s * B + 2mu * A2s / inner_radius^2 * B -
        2 * (lambda + mu) * (1 - theta0) * A1s +
        2mu * (1 - theta0) * A2s / inner_radius^2 - lambda * (1 - theta0) * B +
        (1 - theta0 / 2) * B

    return V0s * pd
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

solversse =
    solver_shell_strain_energy.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0s,
    )
directsse =
    direct_shell_strain_energy.(
        inner_radius,
        outer_radius,
        lambda,
        mu,
        theta0,
        V0s,
    )

err = maximum_error(solversse, directsse)
@test err < 10eps()

solverpotential =
    solver_shell_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0s)
directpotential =
    direct_shell_potential.(inner_radius, outer_radius, lambda, mu, theta0, V0s)

err = maximum_error(solverpotential,directpotential)
