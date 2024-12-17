using FixedPointNumbers,CairoMakie
include("../numerical_scheme.jl")
include("../numerical_scheme3.jl")

#=パラメータの設定
a = 1.0
b = 1.0
c = 1.0
d = 1.0
=#

#=初期値の設定
x = 1.0
y = 1.25
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

function f_float64_rescale_1(x::Float64,y::Float64,t::Float64,scaleratio::Float64)
    r = scaleratio
    a = 1.0
    b = 1.0
    return a*x*r - b*x*y*r^2 
end
function f_float64_rescale_2(x::Float64,y::Float64,t::Float64,scaleratio::Float64)
    r = scaleratio
    c = 1.0
    d = 1.0
    return c*x*y*r^2 - d*y*r
end

function f_BigFloat_1(x::BigFloat,y::BigFloat,t::BigFloat)
    a = BigFloat("1.25")
    b = BigFloat("1.0")
    return a*x - b*x*y 
end
function f_BigFloat_2(x::BigFloat,y::BigFloat,t::BigFloat)
    c = BigFloat("1.0")
    d = BigFloat("1.25")
    return c*x*y - d*y 
end

function f_Q2f61_1(x::Q2f61,y::Q2f61,t::Q2f61)
    a = Q2f61(1.25)
    b = Q2f61(1.0)
    return a*x - b*x*y
end
function f_Q2f61_2(x::Q2f61,y::Q2f61,t::Q2f61)
    c = Q2f61(1.0)
    d = Q2f61(1.25)
    return c*x*y - d*y
end

#Euler法
#Δt = 2^(-7)
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Euler Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Euler Float64 Δt = 2^(-7)")
x = euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-7),Float64,Int(2^(14)))[1]
y = euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-7),Float64,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig

minimum(x)
maximum(x)
minimum(y)
maximum(y)

#rescale_ver
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Euler_rescale Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Euler Float64 Δt = 2^(-7)")
x = euler_rescale_2(f_float64_rescale_1,f_float64_rescale_2,0.5,0.25,2^(-7),1.0,Float64,Int(2^(14)))[1]
y = euler_rescale_2(f_float64_rescale_1,f_float64_rescale_2,0.5,0.25,2^(-7),1.0,Float64,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Euler_rescale_3 Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Euler Float64 Δt = 2^(-7)")
x = euler_rescale_2(f_float64_rescale_1,f_float64_rescale_2,0.5,0.25,2^(-7),3.0,Float64,Int(2^(14)))[1]
y = euler_rescale_2(f_float64_rescale_1,f_float64_rescale_2,0.5,0.25,2^(-7),3.0,Float64,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig

fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra Euler Q2f61 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra Euler Q2f61 Δt = 2^(-7)")
x = euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-7)),Q2f61,Int(2^(14)))[1]
y = euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-7)),Q2f61,Int(2^(14)))[2]
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig

#保存量の計算
function H_float64(x,y)
    a = 1.25
    b = 1.0
    c = 1.0
    d = 1.25
    return c*x + b*y - d*(log2(x)/log2(exp(1))) - a*(log2(y)/log2(exp(1)))
end

function H_BigFloat(x,y)
    a = BigFloat("1.25")
    b = BigFloat("1.0")
    c = BigFloat("1.0")
    d = BigFloat("1.25")
    return c*x + b*y - d*(log2(x)/log2(exp(BigFloat("1.0")))) - a*(log2(y)/log2(exp(BigFloat("1.0"))))
end

function H_Q2f61(x,y)
    a = Q2f61(1.25)
    b = Q2f61(1.0)
    c = Q2f61(1.0)
    d = Q2f61(1.25)
    return c*x + b*y - d*(log2(x)/log2(exp(Q2f61(1.0)))) - a*(log2(y)/log2(exp(Q2f61(1.0))))
end

#Euler法での誤差の比較
#数値解
fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-10,1e2)) ,title = "Lotka-Volterra Euler Numerical Solution Error Δt = 2^(-7)")
ernorm1 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-7),Float64,Int(2^(14)))[1],euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-7),Float64,Int(2^(14)))[2])
#ernorm2 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[2])
ernorm3 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-7)),Q2f61,Int(2^(14)))[1],euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-7)),Q2f61,Int(2^(14)))[2])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax,t,ernorm1,color = :red)
#l2 = lines!(ax,t,ernorm2,color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,ernorm3,color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,#=l2,=#l3],["Float64",#="Q11f52",=#"Q2f61"],position = :rb)
fig

fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "error", yscale = log10, limits = (nothing,(1e-10,1e2)) ,title = "Lotka-Volterra Euler Numerical Solution Error Δt = 2^(-13)")
ernorm1 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[2],euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-13),Float64,Int(2^(20)))[1],euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-13),Float64,Int(2^(20)))[2])
#ernorm2 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-7)),BigFloat,Int(2^(14)))[2],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[2])
ernorm3 = err_norm_2d(euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[1],euler_2(f_BigFloat_1,f_BigFloat_2,BigFloat("1.0"),BigFloat("1.25"),BigFloat(2^(-13)),BigFloat,Int(2^(20)))[2],euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-13)),Q2f61,Int(2^(20)))[1],euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-13)),Q2f61,Int(2^(20)))[2])
t = [i for i in 0:2^(-13):2^7]
l1 = lines!(ax,t,ernorm1,color = :red)
#l2 = lines!(ax,t,ernorm2,color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,ernorm3,color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,#=l2,=#l3],["Float64",#="Q11f52",=#"Q2f61"],position = :rb)
fig

fig = Figure(size =(600, 400))
ax = Axis(fig[1,1],xlabel = "t", ylabel = "H error", yscale = log10, limits = (nothing,(1e-7,1e1)) ,title = "Lotka-Volterra Euler Numerical Solution Error Δt = 2^(-7)")
h_err1 = H_BigFloat(BigFloat("1.0"),BigFloat("1.25")) .- H_float64.(euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-7),Float64,Int(2^(14)))[1],euler_2(f_float64_1,f_float64_2,1.0,1.25,2^(-7),Float64,Int(2^(14)))[2])
#h_err2 = H_BigFloat(BigFloat("1.0"),BigFloat("1.25")) .- H_Q11f52.(euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[1],euler_2(f_Q11f52_1,f_Q11f52_2,Q11f52(1.0),Q11f52(1.25),Q11f52(2^(-7)),Q11f52,Int(2^(14)))[2])
h_err3 = H_BigFloat(BigFloat("1.0"),BigFloat("1.25")) .- H_Q2f61.(euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-7)),Q2f61,Int(2^(14)))[1],euler_2(f_Q2f61_1,f_Q2f61_2,Q2f61(1.0),Q2f61(1.25),Q2f61(2^(-7)),Q2f61,Int(2^(14)))[2])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax,t,abs.(h_err1),color = :red)
#l2 = lines!(ax,t,abs.(h_err2),color = cgrad(:Accent_4)[1])
l3 = lines!(ax,t,abs.(h_err3),color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,#=l2,=#l3],["Float64",#="Q11f52",=#"Q2f61"],position = :rb)
fig