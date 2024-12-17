using FixedPointNumbers,DifferentialEquations,CairoMakie
using Makie.Colors
include("../numerical_scheme.jl")
setprecision(BigFloat,256)
#レスラー方程式
#=a = 0.34375
b = 0.390625
c = 4.5=#
#Float64
function f_fl64_1(x::Float64,y::Float64,z::Float64,t::Float64) 
    return -y - z
end
function f_fl64_2(x::Float64,y::Float64,z::Float64,t::Float64)
    a = 0.34375
    b = 0.390625
    c = 4.5
    return x + a*y
end
function f_fl64_3(x::Float64,y::Float64,z::Float64,t::Float64)
    a = 0.34375
    b = 0.390625
    c = 4.5
    return b + x*z - c*z
end

#BigFloat
function f_BigFloat_1(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    a = BigFloat("0.34375")
    b = BigFloat("0.390625")
    c = BigFloat("4.5")
    return -y - z
end
function f_BigFloat_2(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    a = BigFloat("0.34375")
    b = BigFloat("0.390625")
    c = BigFloat("4.5")
    return x + a*y
end
function f_BigFloat_3(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    a = BigFloat("0.34375")
    b = BigFloat("0.390625")
    c = BigFloat("4.5")
    return b + x*z - c*z
end

#Q8f55
function f_Q8f55_1(x::Q8f55,y::Q8f55,z::Q8f55,t::Q8f55)
    a = Q8f55(0.34375)
    b = Q8f55(0.390625)
    c = Q8f55(4.5)
    return -y - z
end
function f_Q8f55_2(x::Q8f55,y::Q8f55,z::Q8f55,t::Q8f55)
    a = Q8f55(0.34375)
    b = Q8f55(0.390625)
    c = Q8f55(4.5)
    return x + a*y
end
function f_Q8f55_3(x::Q8f55,y::Q8f55,z::Q8f55,t::Q8f55)
    a = Q8f55(0.34375)
    b = Q8f55(0.390625)
    c = Q8f55(4.5)
    return b + x*z - c*z
end

#Q11f52
function f_Q11f52_1(x::Q11f52,y::Q11f52,z::Q11f52,t::Q11f52)
    a = Q11f52(0.34375)
    b = Q11f52(0.390625)
    c = Q11f52(4.5)
    return -y - z
end
function f_Q11f52_2(x::Q11f52,y::Q11f52,z::Q11f52,t::Q11f52)
    a = Q11f52(0.34375)
    b = Q11f52(0.390625)
    c = Q11f52(4.5)
    return x + a*y
end
function f_Q11f52_3(x::Q11f52,y::Q11f52,z::Q11f52,t::Q11f52)
    a = Q11f52(0.34375)
    b = Q11f52(0.390625)
    c = Q11f52(4.5)
    return b + x*z - c*z
end


fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Δt = 0.01")
x = rossler_sq[1]
y = rossler_sq[2]
z = rossler_sq[3]
lines!(ax, x, y, z)  
scatter!(ax, x[1], y[1], z[1])


#==Euler法==#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler BigFloat Δt = 2^(-7)")
x = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1]
y = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2]
z = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3]
lines!(ax1, x, y, z) 
t = [i for i in 0:2^(-7):2^7] 
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler Float64 Δt = 2^(-7)")
x = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[1]
y = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[2]
z = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_float64_euler.pdf", fig) 

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Q8f55 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler Q8f55 Δt = 2^(-7)")
x = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[1]
y = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[2]
z = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_Q8f55_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Q11f52 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler Q11f52 Δt = 2^(-7)")
x = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[1]
y = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[2]   
z = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_Q11f52_euler.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler BigFloat Δt = 2^(-13)")
x = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1]
y = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2]
z = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])  
fig
save("rossler/fig_rossler2/rossler_2^(-13)_BigFloat_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Float64 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler Float64 Δt = 2^(-13)")
x = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[1]
y = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[2]
z = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_float64_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Q8f55 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler Q8f55 Δt = 2^(-13)")
x = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[1]
y = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[2]
z = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_Q8f55_euler.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Q11f52 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Euler Q11f52 Δt = 2^(-13)")
x = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[1]
y = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[2]
z = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_Q11f52_euler.pdf", fig)

#刻み幅　Δt = 1e-6 *BigFloatは時間がかるので分けて行う，まずはfloat64から
fig = Figure(size =(600, 400)) #まだやっていない
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler BigFloat Δt = 1e-6")
x = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),1e-6,BigFloat,Int(1e8))[1]
y = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),1e-6,BigFloat,Int(1e8))[2]
z = euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),1e-6,BigFloat,Int(1e8))[3]
lines!(ax, x, y, z)  
save("rossler/fig_rossler2/rossler_1e-6_BigFloat_euler.pdf", fig)

fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Float64 Δt = 1e-6")
x = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,1e-6,Float64,Int(1e8))[1]
y = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,1e-6,Float64,Int(1e8))[2]
z = euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,1e-6,Float64,Int(1e8))[3]
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_float64_euler.pdf", fig)

fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Q8f55 Δt = 1e-6")
x = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),1e-6,Q8f55,Int(1e8))[1]
y = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),1e-6,Q8f55,Int(1e8))[2]
z = euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),1e-6,Q8f55,Int(1e8))[3]
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_Q8f55_euler.pdf", fig)

fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Euler Q11f52 Δt = 1e-6")
x = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),1e-6,Q11f52,Int(1e8))[1]
y = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),1e-6,Q11f52,Int(1e8))[2]
z = euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),1e-6,Q11f52,Int(1e8))[3]
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_Q11f52_euler.pdf", fig)

#==ルンゲクッタ法==#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta BigFloat Δt = 2^(-7)")
x = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1]
y = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2]
z = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3]
lines!(ax1, x, y, z)  
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta Float64 Δt = 2^(-7)")
x = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[1]
y = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[2]
z = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_float64_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Q8f55 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta Q8f55 Δt = 2^(-7)")
x = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[1]
y = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[2] 
z = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_Q8f55_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Q11f52 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta Q11f52 Δt = 2^(-7)")
x = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[1]
y = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[2]
z = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_Q11f52_rk4.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta BigFloat Δt = 2^(-13)")
x = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1];
y = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2];
z = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3];
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_BigFloat_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Float64 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta Float64 Δt = 2^(-13)")
x = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[1]
y = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[2]
z = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_float64_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Q8f55 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta Q8f55 Δt = 2^(-13)")
x = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[1] 
y = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[2]
z = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_Q8f55_rk4.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Q11f52 Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Runge-Kutta Q11f52 Δt = 2^(-13)")
x = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[1]
y = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[2]
z = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[3]
lines!(ax1, x, y, z)
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], ["x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_Q11f52_rk4.pdf", fig)

#刻み幅　Δt = 1e-6 *BigFloatは時間がかるので分けて行う，まずはfloat64から
fig = Figure(size =(600, 400)) #まだやっていない
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta BigFloat Δt = 1e-6")   
x = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat(4.0),BigFloat(4.0),BigFloat(0.0),1e-6,BigFloat,Int(1e8))[1];
y = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat(4.0),BigFloat(4.0),BigFloat(0.0),1e-6,BigFloat,Int(1e8))[2];
z = rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat(4.0),BigFloat(4.0),BigFloat(0.0),1e-6,BigFloat,Int(1e8))[3];
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_BigFloat_rk4.pdf", fig)

fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Float64 Δt = 1e-6")
x = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,1e-6,Float64,Int(1e8))[1]; 
y = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,1e-6,Float64,Int(1e8))[2];
z = rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,1e-6,Float64,Int(1e8))[3];
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_float64_rk4.pdf", fig)

fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Q8f55 Δt = 1e-6")
x = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),1e-6,Q8f55,Int(1e8))[1];
y = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),1e-6,Q8f55,Int(1e8))[2];
z = rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),1e-6,Q8f55,Int(1e8))[3];
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_Q8f55_rk4.pdf", fig)

fig = Figure(size =(600, 400))
ax = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Runge-Kutta Q11f52 Δt = 1e-6")
x = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),1e-6,Q11f52,Int(1e8))[1];    
y = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),1e-6,Q11f52,Int(1e8))[2];
z = rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),1e-6,Q11f52,Int(1e8))[3];
lines!(ax, x, y, z)
save("rossler/fig_rossler2/rossler_1e-6_Q11f52_rk4.pdf", fig)

#==パッケージを使って厳密解を求める==#
function rossler!(du,u,p,t)
    a,b,c = p
    du[1] = - u[2] - u[3]
    du[2] = u[1] + a*u[2]
    du[3] = b + u[3]*(u[1]-c)
end

tspsan = (BigFloat("0.0"),BigFloat(2^7))
u0 = [BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0")]
p = [BigFloat("0.34375"),BigFloat("0.390625"),BigFloat("4.5")]

prob1 = ODEProblem(rossler!,u0,tspsan,p)
sol1 = solve(prob1,Vern9(),reltol=1e-30,abstol=1e-20,saveat = 2^(-7))
sol2 = solve(prob1,Vern9(),reltol=1e-30,abstol=1e-20,saveat = 2^(-13))

fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Exact Solution BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Exact Solution BigFloat Δt = 2^(-7)")
x = sol1[1,:]
y = sol1[2,:]
z = sol1[3,:]
lines!(ax1, x, y, z)
t = sol1.t
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], [ "x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-7)_BigFloat_exact.pdf", fig)


fig = Figure(size =(800, 400))
ax1 = Axis3(fig[1,1],xlabel = "x",ylabel = "y",zlabel = "z", title = "Rossler Exact Solution BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2],xlabel = "t", title = "Rossler Exact Solution BigFloat Δt = 2^(-13)")
x = sol2[1,:]
y = sol2[2,:]
z = sol2[3,:]
lines!(ax1, x, y, z)
t = sol2.t
l1 = lines!(ax2, t, x)
l2 = lines!(ax2, t, y)
l3 = lines!(ax2, t, z)
axislegend(ax2, [l1, l2, l3], [ "x", "y", "z"])
fig
save("rossler/fig_rossler2/rossler_2^(-13)_BigFloat_exact.pdf", fig)

#刻み幅　Δt = 1e-６はまだ行っていない．

#==Euler法での誤差の比較==#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-16,2^(-7))), title = "Rossler Euler Numerical Solution Error Δt = 2^(-7)")
t_sq = [i for i in 0:2^(-7):2^7]
ernorm1 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[1],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[2],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[3])
ernorm2 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[1],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[2],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[3])
ernorm3 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[3])
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[2])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[1])
axislegend(ax,[l1, l3, l2],["Float64","Q11f52","Q8f55"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-7)_euler_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t",ylabel = "eror", yscale = log10, limits = (nothing,(1e-16,2^(-13))), title = "Rossler Euler Numerical Solution Error Δt = 2^(-13)")
t_sq = [i for i in 0:2^(-13):2^7]
ernorm1 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[1],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[2],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[3])
ernorm2 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[1],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[2],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[3])
ernorm3 = err_norm_3d(euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2],euler_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[3])
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[2])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[1])
axislegend(ax,[l1, l3, l2],["Float64","Q11f52","Q8f55"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-13)_euler_error.pdf", fig)

#刻み幅　Δt = 1e-6 まだやれていない

#==ルンゲクッタ法での誤差の比較==#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-16,2^(-7))), title = "Rossler RK4 Numerical Solution Error Δt = 2^(-7)")
t_sq = [i for i in 0:2^(-7):2^7]
ernorm1 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[1],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[2],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[3])
ernorm2 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[1],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[2],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[3])
ernorm3 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-7),BigFloat,Int(2^14))[3],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[3])
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[2])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[1])
axislegend(ax,[l1, l3, l2],["Float64","Q11f52","Q8f55"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-7)_rk4_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-16,2^(-13))), title = "Rossler RK4 Numerical Solution Error Δt = 2^(-13)")
t_sq = [i for i in 0:2^(-13):2^7]
ernorm1 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[1],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[2],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[3]);
ernorm2 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[1],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[2],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[3]);
ernorm3 = err_norm_3d(rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[1],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[2],rk4_3(f_BigFloat_1,f_BigFloat_2,f_BigFloat_3,BigFloat("4.0"),BigFloat("4.0"),BigFloat("0.0"),2^(-13),BigFloat,Int(2^20))[3],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[3]);
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[2])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[1])
axislegend(ax,[l1, l3, l2],["Float64","Q11f52","Q8f55"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-13)_rk4_error.pdf", fig)

#刻み幅　Δt = 1e-6 までやれていない

#==Euler法と厳密解の比較==#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(2^(-7),1e2)), title = "Rossler Euler Exact Solution Error Δt = 2^(-7)")
t_sq = [i for i in 0:2^(-7):2^7]
ernorm1 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[1],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[2],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[3]);
ernorm2 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[3]);
ernorm3 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[1],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[2],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[3]);
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q8f55","Q11f52"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-7)_euler_exact_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(2^(-13),1e2)), title = "Rossler Euler Exact Solution Error Δt = 2^(-13)")
t_sq = [i for i in 0:2^(-13):2^7]
ernorm1 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[1],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[2],euler_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[3]);
ernorm2 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[1],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[2],euler_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[3]);
ernorm3 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[1],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[2],euler_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[3]);
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q8f55","Q11f52"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-13)_euler_exact_error.pdf", fig)

#刻み幅　Δt = 1e-6 までやれていない

#==RK4法と厳密解の比較==#
#刻み幅　Δt = 2^(-7)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-10,1e2)), title = "Rossler RK4 Exact Solution Error Δt = 2^(-7)")
t_sq = [i for i in 0:2^(-7):2^7]
ernorm1 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[1],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[2],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-7),Float64,Int(2^14))[3]);
ernorm2 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-7),Q11f52,Int(2^14))[3]);
ernorm3 = err_norm_3d(sol1[1,:],sol1[2,:],sol1[3,:],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[1],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[2],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-7),Q8f55,Int(2^14))[3]);
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q8f55","Q11f52"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-7)_rk4_exact_error.pdf", fig)

#刻み幅　Δt = 2^(-13)
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-16,1e2)), title = "Rossler RK4 Exact Solution Error Δt = 2^(-13)")
t_sq = [i for i in 0:2^(-13):2^7]
ernorm1 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[1],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[2],rk4_3(f_fl64_1,f_fl64_2,f_fl64_3,4.0,4.0,0.0,2^(-13),Float64,Int(2^20))[3]);
ernorm2 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[1],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[2],rk4_3(f_Q11f52_1,f_Q11f52_2,f_Q11f52_3,Q11f52(4.0),Q11f52(4.0),Q11f52(0.0),2^(-13),Q11f52,Int(2^20))[3]);
ernorm3 = err_norm_3d(sol2[1,:],sol2[2,:],sol2[3,:],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[1],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[2],rk4_3(f_Q8f55_1,f_Q8f55_2,f_Q8f55_3,Q8f55(4.0),Q8f55(4.0),Q8f55(0.0),2^(-13),Q8f55,Int(2^20))[3]);
l1 = lines!(ax, t_sq, ernorm1, color = :red)
l2 = lines!(ax, t_sq, ernorm2, color = cgrad(:Accent_4)[1])
l3 = lines!(ax, t_sq, ernorm3, color = cgrad(:Accent_4)[2])
axislegend(ax,[l1, l2, l3],["Float64","Q8f55","Q11f52"],position = :lt)
fig
save("rossler/fig_rossler2/rossler_2^(-13)_rk4_exact_error.pdf", fig)