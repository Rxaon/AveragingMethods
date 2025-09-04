using DifferentialEquations, Plots

l = 1; g = 9.81; A = .01; eps = A/l; w = 500; k = g*l/(w^2*A^2); θ0 = 0.1;
t = 10; fps = 60;

y0(t) = sin(t); dy0(t) = cos(t); ddy0(t) = -sin(t);

function f(x, p, t) # x = theta, p, phi, omega
	stheta =
	[
		p[1]*x[2],
		(p[1]*k + ddy0(t))*sin(x[1]),
	]
end

sol = let
	prob = ODEProblem(f, [θ0, 0.0], (0.0, w*t), [eps])
	solve(prob, Vern9(), reltol=1e-8, abstol=1e-8)
end

anim = @animate for ts = 0:w/(fps):w*t
    θ = sol(ts)[1]
    plot([0, l*sin(θ)], [A*y0(ts), A*y0(ts) + l*cos(θ)],
         xlim=[-l,l], ylim=[-l/2,3*l/2], ratio=:equal, legend=false, framestyle=:none)
    scatter!([l*sin(θ)], [A*y0(ts) + l*cos(θ)], ms=20, c=1)
    scatter!([0], [A*y0(ts)], ms=2, c=:black)
end

run(Cmd(`ffmpeg -framerate $fps -pattern_type glob -i '*.png' /home/jorge/movie.mp4`, dir=anim.dir))
