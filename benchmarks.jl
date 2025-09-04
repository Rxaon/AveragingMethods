using DifferentialEquations, LinearAlgebra
import Printf

l = .16; g = 9.81; A = .04; w = 60; k = g*l/(w^2*A^2); gamma = 0.05; θ0 = 0.1;
y0(t) = sin(t); dy0(t) = cos(t); ddy0(t) = -sin(t);

cnt = 0

function f_combined(x, p, t) # x = theta, p, phi, omega
    global cnt; cnt += 1
    stheta = sin(x[1])
	cphi = cos(x[3])
	dphi = p[1]*x[4]/(1+p[1]*y0(t)*cphi)
	sphi = sin(x[3])
	[
		p[1]*x[2],
		(p[1]*k + ddy0(t))*stheta - p[1]*gamma*x[2],
		dphi,
		ddy0(t)*(stheta - sphi) + p[1]*k*stheta - dy0(t)*cphi*dphi - p[1]*gamma*(x[4] + dy0(t)*sphi)
	]
end

function F2(x, p, t)
	int_cond = [x[1] + p[1]*y0(0)*sin(x[1]), x[2] + dy0(0)*sin(x[1]), x[1], x[2]]
	prob_f = ODEProblem(f_combined, int_cond, (0, 2pi), p)
	sol_f = solve(prob_f, Vern9(), reltol=p[2], abstol=p[2])
	prob_b = ODEProblem(f_combined, int_cond, (0, -2pi), p)
	sol_b = solve(prob_b, Vern9(), reltol=p[2], abstol=p[2])
	return (sol_f(2pi)[3:4] - sol_b(-2pi)[3:4])/(4pi)
end

function F4(x, p, t)
	int_cond = [x[1] + p[1]*y0(0)*sin(x[1]), x[2] + dy0(0)*sin(x[1]), x[1], x[2]]
	prob_f = ODEProblem(f_combined, int_cond, (0, 4pi), p)
	sol_f = solve(prob_f, Vern9(), reltol=p[2], abstol=p[2])
	prob_b = ODEProblem(f_combined, int_cond, (0, -4pi), p)
	sol_b = solve(prob_b, Vern9(), reltol=p[2], abstol=p[2])
	return (-sol_f(4pi)[3:4] + 8*(sol_f(2pi)[3:4] - sol_b(-2pi)[3:4]) + sol_b(-4pi)[3:4])/(24pi)
end

function F6(x, p, t)
	int_cond = [x[1] + p[1]*y0(0)*sin(x[1]), x[2] + dy0(0)*sin(x[1]), x[1], x[2]]
	prob_f = ODEProblem(f_combined, int_cond, (0, 6pi), p)
	sol_f = solve(prob_f, Vern9(), reltol=p[2], abstol=p[2])
	prob_b = ODEProblem(f_combined, int_cond, (0, -6pi), p)
	sol_b = solve(prob_b, Vern9(), reltol=p[2], abstol=p[2])
	return (sol_f(6pi)[3:4] - sol_b(-6pi)[3:4] - 9*(sol_f(4pi)[3:4] - sol_b(-4pi)[3:4]) + 45*(sol_f(2pi)[3:4] - sol_b(-2pi)[3:4]))/(120pi)
end

eps = 1e-4; t = 6pi/sqrt(1/2-k);

# Original comp
prob = ODEProblem(f_combined, [θ0, 0.0, θ0, -sin(θ0)], (0.0, t/eps), [eps])
sol = solve(prob, Feagin12(), dt=2pi/7, adaptive=false, saveat=2pi, maxiters=1e8)
println(cnt, ' ')

for tol in [1e-6, 1e-8, 1e-10, 1e-12]
	global cnt

	# Original
	cnt = 0
	prob = ODEProblem(f_combined, [θ0, 0.0, θ0, -sin(θ0)], (0.0, t/eps), [eps])
	sol_new = solve(prob, Vern9(), reltol=tol, abstol=tol, maxiters=1e8)
	s = Printf.@sprintf("%d %.5e ", cnt, maximum(norm(sol(s)[3:4]-sol_new(s)[3:4]) for s in 0:2pi:t/eps))
	s = replace(s, r"e([+-]?)0*(\d+)" => s"e\1\2", '.' => ',')
	print(s)

	# F2
	cnt = 0
	prob = ODEProblem(F2, [θ0, -sin(θ0)], (0.0, t/eps), [eps, tol])
	sol_new = solve(prob, Vern9(), reltol=tol, abstol=tol)
	s = Printf.@sprintf("%d %.5e ", cnt, maximum(norm(sol(s)[3:4]-sol_new(s)) for s in 0:2pi:t/eps))
	s = replace(s, r"e([+-]?)0*(\d+)" => s"e\1\2", '.' => ',')
	print(s)

	# F4
	cnt = 0
	prob = ODEProblem(F4, [θ0, -sin(θ0)], (0.0, t/eps), [eps, tol])
	sol_new = solve(prob, Vern9(), reltol=tol, abstol=tol)
	s = Printf.@sprintf("%d %.5e ", cnt, maximum(norm(sol(s)[3:4]-sol_new(s)) for s in 0:2pi:t/eps))
	s = replace(s, r"e([+-]?)0*(\d+)" => s"e\1\2", '.' => ',')
	print(s)

	# F6
	cnt = 0
	prob = ODEProblem(F6, [θ0, -sin(θ0)], (0.0, t/eps), [eps, tol])
	sol_new = solve(prob, Vern9(), reltol=tol, abstol=tol)
	s = Printf.@sprintf("%d %.5e ", cnt, maximum(norm(sol(s)[3:4]-sol_new(s)) for s in 0:2pi:t/eps))
	s = replace(s, r"e([+-]?)0*(\d+)" => s"e\1\2", '.' => ',')
	println(s)
end
