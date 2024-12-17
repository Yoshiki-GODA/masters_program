using FixedPointNumbers,CairoMakie,DifferentialEquations
include("../numerical_scheme3.jl")
include("../numerical_scheme.jl")

#=パラメータの設定
a = 1.0
b = 1.0
c = 1.0
d = 1.0
=#

#=初期値の設定
x = 0.5
y = 0.25
=#

#保存量 H = cx + by - d(log2(x)/log2(exp(1))) - a(log2(y)/log2(exp(1)))#
function f_float64_1(x::Float64,y::Float64,t::Float64)
    a = 1.0
    b = 1.0
    return a*x - b*x*y
end
function f_float64_2(x::Float64,y::Float64,t::Float64)
    c = 1.0
    d = 1.0
    return c*x*y - d*y
end

function f_BigFloat_1(x::BigFloat,y::BigFloat,t::BigFloat)
    a = BigFloat("1.0")
    b = BigFloat("1.0")
    return a*x - b*x*y  
end
function f_BigFloat_2(x::BigFloat,y::BigFloat,t::BigFloat)
    c = BigFloat("1.0")
    d = BigFloat("1.0")
    return c*x*y - d*y 
end

function f_Q11f52_1(x::Q11f52,y::Q11f52,t::Q11f52)
    a = Q11f52(1.0)
    b = Q11f52(1.0)
    return a*x - b*x*y
end
function f_Q11f52_2(x::Q11f52,y::Q11f52,t::Q11f52)
    c = Q11f52(1.0)
    d = Q11f52(1.0)
    return c*x*y - d*y
end

function f_Q3f60_1(x::Q3f60,y::Q3f60,t::Q3f60)
    a = Q3f60(1.0)
    b = Q3f60(1.0)
    return a*x - b*x*y
end
function f_Q3f60_2(x::Q3f60,y::Q3f60,t::Q3f60)
    c = Q3f60(1.0)
    d = Q3f60(1.0)
    return c*x*y - d*y
end

#保存量の計算
function H_float64(x,y)
    a = 1.0
    b = 1.0
    c = 1.0
    d = 1.0
    return c*x + b*y - d*(log2(x)/log2(exp(1))) - a*(log2(y)/log2(exp(1)))
end

function H_BigFloat(x,y)
    a = BigFloat("1.0")
    b = BigFloat("1.0")
    c = BigFloat("1.0")
    d = BigFloat("1.0")
    return c*x + b*y - d*(log2(x)/log2(exp(BigFloat("1.0")))) - a*(log2(y)/log2(exp(BigFloat("1.0"))))
end

function H_Q11f52(x,y)
    a = Q11f52(1.0)
    b = Q11f52(1.0)
    c = Q11f52(1.0)
    d = Q11f52(1.0)
    return c*x + b*y - d*(log2(x)/log2(exp(Q11f52(1.0)))) - a*(log2(y)/log2(exp(Q11f52(1.0))))
end

function H_Q3f60(x,y)
    a = Q3f60(1.0)
    b = Q3f60(1.0)
    c = Q3f60(1.0)
    d = Q3f60(1.0)
    return c*x + b*y - d*(log2(x)/log2(exp(Q3f60(1.0)))) - a*(log2(y)/log2(exp(Q3f60(1.0))))
end

#stormer_verlet法
#解軌道
#Δt = 2^(-7)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet Float64 Δt = 2^(-7)")
x = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[1]
y = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_sv.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet BigFloat Δt = 2^(-7)")
x = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-7),BigFloat,Int(2^(14)))[1]
y = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-7),BigFloat,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_BigFloat_sv.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet Q11f52 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet Q11f52 Δt = 2^(-7)")
x = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-7),Q11f52,Int(2^(14)))[1]
y = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-7),Q11f52,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q11f52_sv.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet Q3f60 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet Q3f60 Δt = 2^(-7)")
x = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-7),Q3f60,Int(2^(14)))[1]
y = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-7),Q3f60,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q3f60_sv.pdf", fig)

#Δt = 2^(-13)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet Float64 Δt = 2^(-13)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet Float64 Δt = 2^(-13)")
x = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[1]
y = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[2]
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_Float64_sv.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet BigFloat Δt = 2^(-13)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet BigFloat Δt = 2^(-13)")
x = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-13),BigFloat,Int(2^(20)))[1]
y = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-13),BigFloat,Int(2^(20)))[2]
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_BigFloat_sv.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet Q11f52 Δt = 2^(-13)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet Q11f52 Δt = 2^(-13)")
x = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-13),Q11f52,Int(2^(20)))[1]
y = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-13),Q11f52,Int(2^(20)))[2]
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_Q11f52_sv.pdf", fig)

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Stormer Verlet Q3f60 Δt = 2^(-13)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Stormer Verlet Q3f60 Δt = 2^(-13)")
x = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-13),Q3f60,Int(2^(20)))[1]
y = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-13),Q3f60,Int(2^(20)))[2]
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_Q3f60_sv.pdf", fig)

#保存量
H_BigFloat(BigFloat("0.5"),BigFloat("0.25"))
#Δt = 2^(-7)
fig = Figure(size =(600,400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet Float64 Δt = 2^(-7)")
x = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[1]
y = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[2]
h = H_float64.(x,y)
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_sv_h.pdf", fig)

fig = Figure(size =(600,400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet BigFloat Δt = 2^(-7)")
x = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-7),BigFloat,Int(2^(14)))[1]
y = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-7),BigFloat,Int(2^(14)))[2]
h = H_BigFloat.(x,y)
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_BigFloat_sv_h.pdf", fig)

fig = Figure(size =(600,400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet Q11f52 Δt = 2^(-7)")
x = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-7),Q11f52,Int(2^(14)))[1]
y = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-7),Q11f52,Int(2^(14)))[2]
h = H_Q11f52.(x,y)
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q11f52_sv_h.pdf", fig)

fig = Figure(size =(600,400))    
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet Q3f60 Δt = 2^(-7)")
x = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-7),Q3f60,Int(2^(14)))[1]
y = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-7),Q3f60,Int(2^(14)))[2]
h = H_Q3f60.(x,y)
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q3f60_sv_h.pdf", fig)

#Δt = 2^(-13)
fig = Figure(size =(600,400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet Float64 Δt = 2^(-13)")
x = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[1]
y = stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[2]
h = H_float64.(x,y)
t = [i for i in 0:2^(-13):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_float64_sv_h.pdf", fig)

fig = Figure(size =(600,400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet BigFloat Δt = 2^(-13)")
x = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-13),BigFloat,Int(2^(20)))[1]
y = stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),2^(-13),BigFloat,Int(2^(20)))[2]
h = H_BigFloat.(x,y)
t = [i for i in 0:2^(-13):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_BigFloat_sv_h.pdf", fig)

fig = Figure(size =(600,400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet Q11f52 Δt = 2^(-13)")
x = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-13),Q11f52,Int(2^(20)))[1]
y = stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),2^(-13),Q11f52,Int(2^(20)))[2]
h = H_Q11f52.(x,y)
t = [i for i in 0:2^(-13):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_Q11f52_sv_h.pdf", fig)

fig = Figure(size =(600,400))    
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H", title = "Lotka-Volterra Stormer Verlet Q3f60 Δt = 2^(-13)")
x = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-13),Q3f60,Int(2^(20)))[1]
y = stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),2^(-13),Q3f60,Int(2^(20)))[2]
h = H_Q3f60.(x,y)
t = [i for i in 0:2^(-13):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_Q3f60_sv_h.pdf", fig)

#ソルバーを使って解く
function lotka_volterra!(du,u,p,t)
    a,b,c,d = p
    du[1] = a*u[1] - b*u[1]*u[2]
    du[2] = c*u[1]*u[2] - d*u[2]
end

u0 = [BigFloat("0.5"),BigFloat("0.25")]
p = [BigFloat("1.0"),BigFloat("1.0"),BigFloat("1.0"),BigFloat("1.0")]
tspan = (BigFloat("0.0"),BigFloat(2^7))
prob = ODEProblem(lotka_volterra!,u0,tspan,p)
sol1 = solve(prob,Vern9(),reltol=1e-30,abstol=1e-30,saveat = BigFloat(2^(-7)))
sol2 = solve(prob,Vern9(),reltol=1e-30,abstol=1e-30,saveat = BigFloat(2^(-13)))

#stormer_verlet法での誤差の比較
fig = Figure(size =(600,400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error", yscale = log10,limits = (nothing,(1e-20,1e-11)),title = "Lotka-Volterra Stormer Verlet Numerical Solution Error Δt = 2^(-7)")
ernorm1 = err_norm_2d(stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[1],stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[2])
ernorm2 = err_norm_2d(stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[1],stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[2])
ernorm3 = err_norm_2d(stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-7)),Q3f60,Int(2^(14)))[1],stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-7)),Q3f60,Int(2^(14)))[2])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax,t,ernorm1,color = :red)
l2 = lines!(ax,t,ernorm2,color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,ernorm3,color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,l2,l3],["Float64","Q11f52","Q3f60"],position = :rb)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_sv_error.pdf", fig)

fig = Figure(size =(600,400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error", yscale = log10,limits = (nothing,(1e-19,1e-10)),title = "Lotka-Volterra Stormer Verlet Numerical Solution Error Δt = 2^(-13)")
ernorm1 = err_norm_2d(stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[1],stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[2],stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[1],stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[2])
ernorm2 = err_norm_2d(stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[1],stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[2],stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-13)),Q11f52,Int(2^(20)))[1],stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-13)),Q11f52,Int(2^(20)))[2])
ernorm3 = err_norm_2d(stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[1],stormer_verlet(f_BigFloat_1,f_BigFloat_2,BigFloat("0.5"),BigFloat("0.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[2],stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-13)),Q3f60,Int(2^(20)))[1],stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-13)),Q3f60,Int(2^(20)))[2])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax,t,ernorm1,color = :red)
l2 = lines!(ax,t,ernorm2,color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,ernorm3,color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,l2,l3],["Float64","Q11f52","Q3f60"],position = :rb)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_sv_error.pdf", fig)

#保存量
#Δt = 2^(-7)
fig = Figure(size =(600,400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error", yscale = log10,limits = (nothing,(1e-10,1e-1)),title = "Lotka-Volterra Stormer Verlet Numerical Solution Error Δt = 2^(-7)")
h_err1 = H_BigFloat(BigFloat("0.5"),BigFloat("0.25")) .- H_float64.(stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[1],stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-7),Float64,Int(2^(14)))[2])
h_err2 = H_BigFloat(BigFloat("0.5"),BigFloat("0.25")) .- H_Q11f52.(stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[1],stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[2])
h_err3 = H_BigFloat(BigFloat("0.5"),BigFloat("0.25")) .- H_Q3f60.(stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-7)),Q3f60,Int(2^(14)))[1],stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-7)),Q3f60,Int(2^(14)))[2])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax,t,abs.(h_err1),color = :red)
l2 = lines!(ax,t,abs.(h_err2),color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,abs.(h_err3),color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,l2,l3],["Float64","Q11f52","Q3f60"],position = :rb)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_sv_h_error.pdf", fig)

#Δt = 2^(-13)
fig = Figure(size =(600,400))
ax = Axis(fig[1,1], xlabel = "t", ylabel = "error", yscale = log10,limits = (nothing,(1e-13,1e-4)),title = "Lotka-Volterra Stormer Verlet Numerical Solution Error Δt = 2^(-13)")   
h_err1 = H_BigFloat(BigFloat("0.5"),BigFloat("0.25")) .- H_float64.(stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[1],stormer_verlet(f_float64_1,f_float64_2,0.5,0.25,2^(-13),Float64,Int(2^(20)))[2])
h_err2 = H_BigFloat(BigFloat("0.5"),BigFloat("0.25")) .- H_Q11f52.(stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-13)),Q11f52,Int(2^(20)))[1],stormer_verlet(f_Q11f52_1,f_Q11f52_2,Q11f52(0.5),Q11f52(0.25),Q11f52(2^(-13)),Q11f52,Int(2^(20)))[2])
h_err3 = H_BigFloat(BigFloat("0.5"),BigFloat("0.25")) .- H_Q3f60.(stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-13)),Q3f60,Int(2^(20)))[1],stormer_verlet(f_Q3f60_1,f_Q3f60_2,Q3f60(0.5),Q3f60(0.25),Q3f60(2^(-13)),Q3f60,Int(2^(20)))[2])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax,t,abs.(h_err1),color = :red)
l2 = lines!(ax,t,abs.(h_err2),color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,abs.(h_err3),color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,l2,l3],["Float64","Q11f52","Q3f60"],position = :rb)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-13)_sv_h_error.pdf", fig)