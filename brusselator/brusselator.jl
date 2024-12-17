using FixedPointNumbers, DifferentialEquations, CairoMakie,Makie.Colors
include("../numerical_scheme.jl")

#Brusselator方程式
#=  a = 1
    b = 3=#
#Float64
function f_fl64_1(x::Float64,y::Float64,t::Float64)
    a = 1.0
    b = 3.0
    return a + x^2*y - b*x - x
end
function f_fl64_2(x::Float64,y::Float64,t::Float64)
    a = 1.0
    b = 3.0
    return -(x^2)*y + b*x
end

#BigFloat
function f_BigFloat_1(x::BigFloat,y::BigFloat,t::BigFloat)
    a = BigFloat("1.0")
    b = BigFloat("3.0")
    return a + x^2*y - b*x - x
end
function f_BigFloat_2(x::BigFloat,y::BigFloat,t::BigFloat)
    a = BigFloat("1.0")
    b = BigFloat(3.0)
    return -(x^2)*y + b*x
end

#Q11f52
function f_Q11f52_1(x::Q11f52,y::Q11f52,t::Q11f52)
    a = Q11f52(1.0)
    b = Q11f52(3.0)
    return a + x^2*y - b*x - x
end
function f_Q11f52_2(x::Q11f52,y::Q11f52,t::Q11f52)
    a = Q11f52(1.0)
    b = Q11f52(3.0)
    return -(x^2)*y + b*x
end

function f_Q7f56_1(x::Q7f56,y::Q7f56,t::Q7f56)
    a = Q7f56(1.0)
    b = Q7f56(3.0)
    return a + x^2*y - b*x - x
end
function f_Q7f56_2(x::Q7f56,y::Q7f56,t::Q7f56)
    a = Q7f56(1.0)
    b = Q7f56(3.0)
    return -(x^2)*y + b*x
end

function f_Q6f57_1(x::Q6f57,y::Q6f57,t::Q6f57)
    a = Q6f57(1.0)
    b = Q6f57(3.0)
    return a + x^2*y - b*x - x   
end
function f_Q6f57_2(x::Q6f57,y::Q6f57,t::Q6f57)
    a = Q6f57(1.0)
    b = Q6f57(3.0)
    return -(x^2)*y + b*x
end

function f_Q5f58_1(x::Q5f58,y::Q5f58,t::Q5f58)
    a = Q5f58(1.0)
    b = Q5f58(3.0)
    return a + x^2*y - b*x - x
end
function f_Q5f58_2(x::Q5f58,y::Q5f58,t::Q5f58)
    a = Q5f58(1.0)
    b = Q5f58(3.0)
    return -(x^2)*y + b*x
end

function f_Q5f58_1(x::Q5f58,y::Q5f58,t::Q5f58)
    a = Q5f58(1.0)
    b = Q5f58(3.0)
    return a + x^2*y - b*x - x
end
function f_Q5f58_2(x::Q5f58,y::Q5f58,t::Q5f58)
    a = Q5f58(1.0)
    b = Q5f58(3.0)
    return -(x^2)*y + b*x
end

function f_Q4f59_1(x::Q4f59,y::Q4f59,t::Q4f59)
    a = Q4f59(1.0)
    b = Q4f59(3.0)
    return a + x^2*y - b*x - x   
end
function f_Q4f59_2(x::Q4f59,y::Q4f59,t::Q4f59)
    a = Q4f59(1.0)
    b = Q4f59(3.0)
    return -(x^2)*y + b*x
end

function f_Q3f60_1(x::Q3f60,y::Q3f60,t::Q3f60)
    a = Q3f60(1.0)
    b = Q3f60(3.0)
    return a + x^2*y - b*x - x
end
function f_Q3f60_2(x::Q3f60,y::Q3f60,t::Q3f60)
    a = Q3f60(1.0)
    b = Q3f60(3.0)
    return -(x^2)*y + b*x
end

#=Euler法=#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1],xlabel = "x", ylabel = "y", title = "Brusselator Euler BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler BigFloat Δt = 2^(-7)")
x = euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1]
y = euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2]
t = [i for i in 0:2^(-7):2^5]
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1],xlabel = "x", ylabel = "y", title = "Brusselator Euler Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Float64 Δt = 2^(-7)")
x = euler_2(f_fl64_1,f_fl64_2,Float64(1.0),Float64(1.0),Float64(2^(-7)),Float64,Int(2^12))[1];
y = euler_2(f_fl64_1,f_fl64_2,Float64(1.0),Float64(1.0),Float64(2^(-7)),Float64,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_float64_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1],xlabel = "x", ylabel = "y", title = "Brusselator Euler Q11f52 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Q11f52 Δt = 2^(-7)")
x = euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[1];
y = euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_Q11f52_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1],xlabel = "x", ylabel = "y", title = "Brusselator Euler Q4f59 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Q4f59 Δt = 2^(-7)")
x = euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[1];
y = euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_Q4f59_euler.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1],xlabel = "x", ylabel = "y", title = "Brusselator Euler BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler BigFloat Δt = 2^(-13)")
x = euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1];
y = euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler Float64 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Float64 Δt = 2^(-13)")
x = euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[1];
y = euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_float64_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler Q11f52 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Q11f52 Δt = 2^(-13)")
x = euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[1];
y = euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_Q11f52_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler Q4f59 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Q4f59 Δt = 2^(-13)")
x = euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[1];
y = euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_Q4f59_euler.pdf", fig)

#刻み幅　Δt = 1e-6
fig = Figure(size =(800, 400)) #BigFloatはダメでした．．．
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler BigFloat Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler BigFloat Δt = 1e-6")
x = euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(1e-6),BigFloat,Int(3e7))[1];
y = euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(1e-6),BigFloat,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler Float64 Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Float64 Δt = 1e-6")
x = euler_2(f_fl64_1,f_fl64_2,1.0,1.0,1e-6,Float64,Int(3e7))[1];
y = euler_2(f_fl64_1,f_fl64_2,1.0,1.0,1e-6,Float64,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [elm1, elm2], [ "x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_float64_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler Q11f52 Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Q11f52 Δt = 1e-6")
x = euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(1e-6),Q11f52,Int(3e7))[1];
y = euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(1e-6),Q11f52,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [ elm1, elm2], ["x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_Q11f52_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Euler Q4f59 Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Euler Q4f59 Δt = 1e-6")
x = euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(1e-6),Q4f59,Int(3e7))[1];
y = euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(1e-6),Q4f59,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [elm1, elm2], ["x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_Q4f59_euler.pdf", fig)

#=ルンゲ・クッタ法=#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q11f52 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q11f52 Δt = 2^(-7)")
x = rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1];
y = rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Float64 Δt = 2^(-7)")
x = rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[1];
y = rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_float64_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q11f52 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q11f52 Δt = 2^(-7)")
x = rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[1];
y = rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_Q11f52_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q4f59 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q4f59 Δt = 2^(-7)")
x = rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[1];
y = rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[2];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_Q4f59_rk4.pdf", fig)

#刻み幅：Δt = 2^(-13)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta BigFloat Δt = 2^(-13)")
x = rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1];
y = rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Float64 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Float64 Δt = 2^(-13)")
x = rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[1];  
y = rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_float64_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q11f52 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q11f52 Δt = 2^(-13)")
x = rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[1];
y = rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_Q11f52_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q4f59 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q4f59 Δt = 2^(-13)")
x = rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[1];
y = rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[2];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_Q4f59_rk4.pdf", fig)

#刻み幅：Δt = 1e-6
fig = Figure(size =(800, 400)) #まだやていない．．．
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta BigFloat Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta BigFloat Δt = 1e-6")
x = rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(1e-6),BigFloat,Int(3e7))[1];
y = rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(1e-6),BigFloat,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [elm1, elm2], ["x", "y"])

save("brusselator/fig_brusselator2/brusselator_1e-6_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Float64 Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Float64 Δt = 1e-6")
x = rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,1e-6,Float64,Int(3e7))[1];
y = rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,1e-6,Float64,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [elm1, elm2], ["x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_float64_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q11f52 Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q11f52 Δt = 1e-6")
x = rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(1e-6),Q11f52,Int(3e7))[1];
y = rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(1e-6),Q11f52,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [elm1, elm2], ["x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_Q11f52_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Runge-Kutta Q4f59 Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Runge-Kutta Q4f59 Δt = 1e-6")
x = rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(1e-6),Q4f59,Int(3e7))[1];
y = rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(1e-6),Q4f59,Int(3e7))[2];
t = [i for i in 0:1e-6:2^5];
l1 = scatter!(ax1, x, y, markersize = 0.5)
l2 = scatter!(ax2, t, x, markersize = 0.5)
l3 = scatter!(ax2, t, y, markersize = 0.5)
elm1 = LineElement(color = 0,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
elm2 = LineElement(color = 1,linestyle = nothing,colorrange = (0, 10),colormap = :tab10)
Legend(fig[1,3], [elm1, elm2], ["x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_Q4f59_rk4.pdf", fig)

#=パッケージで厳密解を求める=#
function Brusselator!(du,u,p,t)
    a,b = p
    du[1] = a + (u[1]^2)*u[2] - b*u[1] - u[1]
    du[2] = b*u[1] - (u[1]^2)*u[2]
end

tspan = (BigFloat("0.0"),BigFloat(2^5))
u0 = [BigFloat("1.0");BigFloat("1.0")]
p = [BigFloat("1.0"),BigFloat("3.0")]

prob1 = ODEProblem(Brusselator!,u0,tspan,p)
sol1 = solve(prob1,Vern9(),reltol=1e-30,abstol=1e-20,saveat = 2^(-7))
sol2 = solve(prob1,Vern9(),reltol=1e-30,abstol=1e-20,saveat = 2^(-13))
#sol3 = solve(prob1,Vern9(),reltol=1e-8,abstol=1e-8,saveat = 1e-6)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Exact Solution BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Exact Solution BigFloat Δt = 2^(-7)")
x = sol1[1,:];
y = sol1[2,:];
t = [i for i in 0:2^(-7):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_BigFloat_exact.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Exact Solution BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Exact Solution BigFloat Δt = 2^(-13)")
x = sol2[1,:];
y = sol2[2,:];
t = [i for i in 0:2^(-13):2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_BigFloat_exact.pdf", fig)

fig = Figure(size =(800, 400)) #まだできていない．．．
ax1 = Axis(fig[1,1], xlabel = "x", ylabel = "y", title = "Brusselator Exact Solution BigFloat Δt = 1e-6")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Brusselator Exact Solution BigFloat Δt = 1e-6")
x = sol3[1,:];
y = sol3[2,:];
t = [i for i in 0:1e-6:2^5];
l1 = lines!(ax1, x, y)
l2 = lines!(ax2, t, x)
l3 = lines!(ax2, t, y)
axislegend(ax2, [ l2, l3], [ "x", "y"])
save("brusselator/fig_brusselator2/brusselator_1e-6_BigFloat_exact.pdf", fig)

#=オイラー法での誤差の比較=#
#刻み幅：Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-23,1e-12)), title = "Brusselator Euler Numerical Solution Error Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^5];
ernorm1 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[1],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[2]);
ernorm2 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[2]);
ernorm3 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[1],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :rb)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_euler_error.pdf", fig)

#刻み幅：Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-21,1e-10)), title = "Brusselator Euler Numerical Solution Error Δt = 2^(-13)")
t = [i for i in 0:2^(-13):2^5];
ernorm1 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[1],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[2]);
ernorm2 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[2]);
ernorm3 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[1],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :rb)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_euler_error.pdf", fig)

#=ルンゲ・クッタ法での誤差の比較=#
#刻み幅：Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-21,1e-12)), title = "Brusselator Runge-Kutta Numerical Solution Error Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^5];
ernorm1 = err_norm_2d(rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1],rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[1],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[2]);
ernorm2 = err_norm_2d(rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1],rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[1],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[2]);
ernorm3 = err_norm_2d(rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[1],rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-7)),BigFloat,Int(2^12))[2],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[1],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :rb)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_rk4_error.pdf", fig)

#刻み幅：Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-21,1e-10)), title = "Brusselator Runge-Kutta Numerical Solution Error Δt = 2^(-13)")
t = [i for i in 0:2^(-13):2^5];
ernorm1 = err_norm_2d(rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1],rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[1],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[2]);
ernorm2 = err_norm_2d(rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1],rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[1],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[2]);
ernorm3 = err_norm_2d(rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[1],rk4_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.0"),BigFloat(2^(-13)),BigFloat,Int(2^18))[2],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[1],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :rb)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_rk4_error.pdf", fig)

#=オイラー法と厳密解の比較=#
#刻み幅：Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-5,1e1)), title = "Brusselator Euler Exact Solution Error Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^5];
ernorm1 = err_norm_2d(sol1[1,:],sol1[2,:],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[1],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[2]);
ernorm2 = err_norm_2d(sol1[1,:],sol1[2,:],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-7)),Q11f52,Int(2^12))[2]);
ernorm3 = err_norm_2d(sol1[1,:],sol1[2,:],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[1],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-7)),Q4f59,Int(2^12))[2]);
l1 = lines!(ax, t, ernorm1, color = :red) 
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :rb)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_euler_exact_error.pdf", fig)

#刻み幅：Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-7,1e-1)), title = "Brusselator Euler Exact Solution Error Δt = 2^(-13)")
t = [i for i in 0:2^(-13):2^5];
ernorm1 = err_norm_2d(sol2[1,:],sol2[2,:],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[1],euler_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[2]);
ernorm2 = err_norm_2d(sol2[1,:],sol2[2,:],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),Q11f52(2^(-13)),Q11f52,Int(2^18))[2]);
ernorm3 = err_norm_2d(sol2[1,:],sol2[2,:],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[1],euler_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),Q4f59(2^(-13)),Q4f59,Int(2^18))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :lt)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_euler_exact_error.pdf", fig)

#=ルンゲ・クッタ法と厳密解の比較=#
#刻み幅：Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-12,1e-5)), title = "Brusselator Runge-Kutta Exact Solution Error Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^5];
ernorm1 = err_norm_2d(sol1[1,:],sol1[2,:],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[1],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-7),Float64,Int(2^12))[2]);
ernorm2 = err_norm_2d(sol1[1,:],sol1[2,:],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),2^(-7),Q11f52,Int(2^12))[1],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),2^(-7),Q11f52,Int(2^12))[2]);
ernorm3 = err_norm_2d(sol1[1,:],sol1[2,:],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),2^(-7),Q4f59,Int(2^12))[1],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),2^(-7),Q4f59,Int(2^12))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :lt)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-7)_rk4_exact_error.pdf", fig)

#刻み幅：Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error",yscale = log10, limits = (nothing,(1e-19,1e-6)), title = "Brusselator Runge-Kutta Exact Solution Error Δt = 2^(-13)")
t = [i for i in 0:2^(-13):2^5];
ernorm1 = err_norm_2d(sol2[1,:],sol2[2,:],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[1],rk4_2(f_fl64_1,f_fl64_2,1.0,1.0,2^(-13),Float64,Int(2^18))[2]);
ernorm2 = err_norm_2d(sol2[1,:],sol2[2,:],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),2^(-13),Q11f52,Int(2^18))[1],rk4_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.0),2^(-13),Q11f52,Int(2^18))[2]);
ernorm3 = err_norm_2d(sol2[1,:],sol2[2,:],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),2^(-13),Q4f59,Int(2^18))[1],rk4_2(f_Q4f59_1,f_Q4f59_2,Q4f59(1.0),Q4f59(1.0),2^(-13),Q4f59,Int(2^18))[2]);
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_5)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_5)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q4f59"],position = :rb)
fig
save("brusselator/fig_brusselator2/brusselator_2^(-13)_rk4_exact_error.pdf", fig)