using CairoMakie, FixedPointNumbers, DifferentialEquations
CairoMakie.activate!(type = "pdf")
CairoMakie.activate!(inline = true)
include("../numerical_scheme.jl")
setprecision(BigFloat,256)
sizeof(BigFloat)*8 - precision(BigFloat) -1
#Lorenz方程式

#=
σ = 10.0
ρ = 28.0
β = 2.6640625
=#

#刻み幅2^(-7)

#float64
function f_float64_1(x::Float64,y::Float64,z::Float64,t::Float64)
    σ = 10.0
    return σ*(y - x)
end
function f_float64_2(x::Float64,y::Float64,z::Float64,t::Float64)
    ρ = 28.0
    return x*(ρ - z) - y
end
function f_float64_3(x::Float64,y::Float64,z::Float64,t::Float64)
    β = 2.6640625
    return x*y - β*z
end

#64ビットの固定小数点数で実行
function f_Q11f52_1(x::Q11f52,y::Q11f52,z::Q11f52,t::Q11f52)
    σ_Q11f52 = Q11f52(10.0)
    return σ_Q11f52*(y - x)
end
function f_Q11f52_2(x::Q11f52,y::Q11f52,z::Q11f52,t::Q11f52)
    ρ_Q11f52 = Q11f52(28.0)
    return x*(ρ_Q11f52 - z) - y
end
function f_Q11f52_3(x::Q11f52,y::Q11f52,z::Q11f52,t::Q11f52)
    β_Q11f52 = Q11f52(2.6640625)
    return x*y - β_Q11f52*z
end

function f_Q9f54_1(x::Q9f54,y::Q9f54,z::Q9f54,t::Q9f54)
    σ_Q9f54 = Q9f54(10.0)
    return σ_Q9f54*(y - x)    
end
function f_Q9f54_2(x::Q9f54,y::Q9f54,z::Q9f54,t::Q9f54)
    ρ_Q9f54 = Q9f54(28.0)
    return x*(ρ_Q9f54 - z) - y
end
function f_Q9f54_3(x::Q9f54,y::Q9f54,z::Q9f54,t::Q9f54)
    β_Q9f54 = Q9f54(2.6640625)
    return x*y - β_Q9f54*z
end

#BigFloatで実行
function f_BigFloat_1(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    σ_BigFloat = BigFloat("10.0")
    return σ_BigFloat*(y - x)
end
function f_BigFloat_2(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    ρ_BigFloat = BigFloat("28.0")
    return x*(ρ_BigFloat - z) - y
end
function f_BigFloat_3(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    β_BigFloat = BigFloat("2.6640625")
    return x*y - β_BigFloat*z
end

#===Euler法===#
#刻み幅　Δt = 2^(-7)
#数値解を図示
fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler BigFloat, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t", ylabel = "x",title = "Lorenz Euler BigFloat, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t", ylabel = "y",title = "Lorenz Euler BigFloat, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t", ylabel = "z",title = "Lorenz Euler BigFloat, Δt = 2^(-7)")
x = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1]
y = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2]
z = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
fig
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
save("lorenz/fig_lorenz/lorenz_2^(-7)_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler float64, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t", ylabel = "x" ,title = "Lorenz Euler float64, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t", ylabel = "y" ,title = "Lorenz Euler float64, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t", ylabel = "z" ,title = "Lorenz Euler float64, Δt = 2^(-7)")
x = euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[1]
y = euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[2]
z = euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_float64_euler.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler Q11f52, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Euler Q11f52, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Euler Q11f52, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Euler Q11f52, Δt = 2^(-7)")
x = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-7)),Q11f52,Int(2^14))[1]
y = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-7)),Q11f52,Int(2^14))[2]
z = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-7)),Q11f52,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_Q11f52_euler.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler Q9f54, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Euler Q9f54, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Euler Q9f54, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Euler Q9f54, Δt = 2^(-7)")
x = euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-7)),Q9f54,Int(2^14))[1]
y = euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-7)),Q9f54,Int(2^14))[2]
z = euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-7)),Q9f54,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_Q9f54_euler.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler BigFloat, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Euler BigFloat, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Euler BigFloat, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Euler BigFloat, Δt = 2^(-13)")
x = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1]
y = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2]
z = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler float64, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Euler float64, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Euler float64, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Euler float64, Δt = 2^(-13)")
x = euler_3(f_float64_1,f_float64_2,f_float64_3,Float64(0.5),Float64(0.5),Float64(0.5),Float64(2^(-13)),Float64,Int(2^20))[1]
y = euler_3(f_float64_1,f_float64_2,f_float64_3,Float64(0.5),Float64(0.5),Float64(0.5),Float64(2^(-13)),Float64,Int(2^20))[2]
z = euler_3(f_float64_1,f_float64_2,f_float64_3,Float64(0.5),Float64(0.5),Float64(0.5),Float64(2^(-13)),Float64,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_float64_euler.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler Q11f52, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Euler Q11f52, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Euler Q11f52, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Euler Q11f52, Δt = 2^(-13)")
x = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-13)),Q11f52,Int(2^20))[1]
y = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-13)),Q11f52,Int(2^20))[2]
z = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-13)),Q11f52,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_Q11f52_euler.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Euler Q9f54, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Euler Q9f54, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Euler Q9f54, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Euler Q9f54, Δt = 2^(-13)")
x = euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-13)),Q9f54,Int(2^20))[1]
y = euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-13)),Q9f54,Int(2^20))[2]
z = euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-13)),Q9f54,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_Q9f54_euler.pdf", fig)

#===ルンゲ・クッタ法===#
#刻み幅　Δt = 2^(-7)
#数値解を図示
fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-7)")
x = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1]
y = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2]
z = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta float64, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta float64, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta float64, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta float64, Δt = 2^(-7)")
x = rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[1]
y = rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[2]
z = rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_float64_rk4.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-7)")
x = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-7)),Q11f52,Int(2^14))[1]
y = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-7)),Q11f52,Int(2^14))[2]
z = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-7)),Q11f52,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_Q11f52_rk4.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-7)")
x = rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-7)),Q9f54,Int(2^14))[1]
y = rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-7)),Q9f54,Int(2^14))[2]
z = rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-7)),Q9f54,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_Q9f54_rk4.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t", ylabel = "x",title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t", ylabel = "y",title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t", ylabel = "z",title = "Lorenz Runge-Kutta BigFloat, Δt = 2^(-13)")
x = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1]
y = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2]
z = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta Float64, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta Float64, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta Float64, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta Float64, Δt = 2^(-13)")
x = rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[1]
y = rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[2]
z = rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)    
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_Float64_rk4.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta Q11f52, Δt = 2^(-13)")
x = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-13)),Q11f52,Int(2^20))[1]
y = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-13)),Q11f52,Int(2^20))[2]
z = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),Q11f52(2^(-13)),Q11f52,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)    
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_Q11f52_rk4.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Runge-Kutta Q9f54, Δt = 2^(-13)")
x = rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-13)),Q9f54,Int(2^20))[1]
y = rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-13)),Q9f54,Int(2^20))[2]
z = rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),Q9f54(2^(-13)),Q9f54,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)    
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_Q9f54_rk4.pdf", fig)

#===微分方程式ソルバーを実行===#
function Lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

tspan = (BigFloat("0.0"),BigFloat("128"))
u0 = [BigFloat("0.5");BigFloat("0.5");BigFloat("0.5")]
p = [BigFloat("10.0"),BigFloat("28.0"),BigFloat("2.6640625")]

prob = ODEProblem(Lorenz!,u0,tspan,p)
sol1 = solve(prob,Vern9(),reltol=1e-30,abstol=1e-20,saveat = BigFloat(2^(-7)))
sol2 = solve(prob,Vern9(),reltol=1e-30,abstol=1e-20,saveat = BigFloat(2^(-13)))

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Exact Solution BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Exact Solution BigFloat Δt = 2^(-7)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Exact Solution BigFloat Δt = 2^(-7)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Exact Solution BigFloat Δt = 2^(-7)")
x = sol1[1,:]
y = sol1[2,:]
z = sol1[3,:]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)    
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_BigFloat_exact.pdf", fig)

fig = Figure(size =(800, 600))
ax1 = Axis3(fig[1, 1],xlabel ="x", ylabel = "y", zlabel = "z", viewmode = :fit, xreversed = true, title = "Lorenz Exact Solution BigFloatΔt = 2^(-13)")
ax2 = Axis(fig[1,2][1,1],xlabel = "t",ylabel = "x", title = "Lorenz Exact Solution BigFloat Δt = 2^(-13)")
ax3 = Axis(fig[1,2][2,1],xlabel = "t",ylabel = "y", title = "Lorenz Exact Solution BigFloat Δt = 2^(-13)")
ax4 = Axis(fig[1,2][3,1],xlabel = "t",ylabel = "z", title = "Lorenz Exact Solution BigFloat Δt = 2^(-13)")
x = sol2[1,:]
y = sol2[2,:]
z = sol2[3,:]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2,t,x)
l2 = lines!(ax3,t,y)    
l3 = lines!(ax4,t,z)
#axislegend(ax2,[l1,l2,l3],["x","y","z"])
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_BigFloat_exact.pdf", fig)

#=Euler法での誤差の比較=#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-15,1e2)), title = "Lorenz Euler Numerical Solution Error Δt = 2^(-7)")
ernorm1 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[1],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[2],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[3])
ernorm2 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[3]) 
ernorm3 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[1],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[2],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[3])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_euler_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-15,1e2)), title = "Lorenz Euler Numerical Solution Error Δt = 2^(-13)")
ernorm1 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[1],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[2],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[3])
ernorm2 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[3]) 
ernorm3 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[1],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[2],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[3])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_euler_error.pdf", fig)

#=ルンゲ・クッタ法での誤差の比較=#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-17,1e2)), title = "Lorenz Runge-Kutta Numerical Solution Error Δt = 2^(-7)")
ernorm1 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[1],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[2],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[3])
ernorm2 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[3]) 
ernorm3 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-7)),BigFloat,Int(2^14))[3],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[1],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[2],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[3])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_rk4_error.pdf", fig)

#刻み幅 Δt = 2^(-13)
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-17,1e2)), title = "Lorenz Runge-Kutta Numerical Solution Error Δt = 2^(-13)")
ernorm1 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[1],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[2],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[3])
ernorm2 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[3]) 
ernorm3 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("0.5"),BigFloat("0.5"),BigFloat("0.5"),BigFloat(2^(-13)),BigFloat,Int(2^20))[3],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[1],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[2],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[3])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_rk4_error.pdf", fig)

#=Euler法と厳密解の誤差=#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-5,1e2)), title = "Lorenz Euler Exact Solution Error Δt = 2^(-7)")
ernorm1 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[1],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[2],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[3])
ernorm2 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[3]) 
ernorm3 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[1],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[2],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[3])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_euler_exact_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-8,1e2)), title = "Lorenz Euler Exact Solution Error Δt = 2^(-13)")
ernorm1 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[1],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[2],euler_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[3])
ernorm2 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[3]) 
ernorm3 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[1],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[2],euler_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[3])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_euler_exact_error.pdf", fig)
#=ルンゲ・クッタ法と厳密解との誤差=#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-5,1e2)), title = "Lorenz Runge-Kutta Exact Solution Error Δt = 2^(-7)")
ernorm1 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[1],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[2],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-7),Float64,Int(2^14))[3])
ernorm2 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-7),Q11f52,Int(2^14))[3]) 
ernorm3 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[1],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[2],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-7),Q9f54,Int(2^14))[3])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-7)_rk4_exaxt_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1],xlabel ="t", ylabel = "error", yscale = log10, limits = (nothing,(1e-8,1e2)), title = "Lorenz Runge-Kutta Exact Solution Error Δt = 2^(-13)")
ernorm1 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[1],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[2],rk4_3(f_float64_1,f_float64_2,f_float64_3,0.5,0.5,0.5,2^(-13),Float64,Int(2^20))[3])
ernorm2 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(0.5),Q11f52(0.5),Q11f52(0.5),2^(-13),Q11f52,Int(2^20))[3]) 
ernorm3 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[1],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[2],rk4_3(f_Q9f54_1,f_Q9f54_2,f_Q9f54_3,Q9f54(0.5),Q9f54(0.5),Q9f54(0.5),2^(-13),Q9f54,Int(2^20))[3])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax, t, ernorm1, color = :red)
l2 = lines!(ax, t, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q11f52","Q9f54"],position = :rb)
fig
save("lorenz/fig_lorenz/lorenz_2^(-13)_rk4_exact_error.pdf", fig)