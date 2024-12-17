using LinearAlgebra,FixedPointNumbers,CairoMakie
include("../numerical_scheme.jl")
include("../numerical_scheme3.jl")
setprecision(BigFloat,256)

#微分方程式
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

function f_Q4f59_1(x::Q4f59,y::Q4f59,t::Q4f59)
    a = Q4f59(1.0)
    b = Q4f59(1.0)
    return a*x - b*x*y
end

function f_Q4f59_2(x::Q4f59,y::Q4f59,t::Q4f59)
    c = Q4f59(1.0)
    d = Q4f59(1.0)
    return c*x*y - d*y
end

#spmの更新式のベクトル
function vec_spm_BigFloat(v_new::Vector{BigFloat},v_old::Vector{BigFloat},s::BigFloat,p::Array{BigFloat,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ]
end

function vec_spm_float64(v_new::Vector{Float64},v_old::Vector{Float64},s::Float64,p::Vector{Float64})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ]
end

function vec_spm_Q11f52(v_new::Vector{Q11f52},v_old::Vector{Q11f52},s::Q11f52,p::Vector{Q11f52})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        Q11f52(x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old))),
        Q11f52(y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old)))
    ]
end

function vec_spm_Q4f59(v_new::Vector{Q4f59},v_old::Vector{Q4f59},s::Q4f59,p::Array{Q4f59,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        Q4f59(x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old))),
        Q4f59(y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old)))
    ] 
end

#更新式の行列
function mat_spm_BigFloat(v_new::Vector{BigFloat},v_old::Vector{BigFloat},s::BigFloat,p::Vector{BigFloat})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        1 s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2)))
        s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2))) 1
    ]
end

function mat_spm_float64(v_new::Vector{Float64},v_old::Vector{Float64},s::Float64,p::Vector{Float64})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        1 s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2)))
        s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2))) 1
    ]
end

function mat_spm_Q11f52(v_new::Vector{Q11f52},v_old::Vector{Q11f52},s::Q11f52,p::Vector{Q11f52})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        1 Q11f52(s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2))))
        Q11f52(s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2)))) 1
    ]
end

function mat_spm_Q4f59(v_new::Vector{Q4f59},v_old::Vector{Q4f59},s::Q4f59,p::Vector{Q4f59})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        1 Q4f59(s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2))))
        Q4f59(s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2)))) 1
    ]
end

#保存量の計算
function H_BigFloat_exp(x,y)
    a = BigFloat("1.0")
    b = BigFloat("1.0")
    c = BigFloat("1.0")
    d = BigFloat("1.0")
    return c*exp(x) + b*exp(y) - d*x - a*y
end

function H_float64_exp(x,y)
    a = 1.0
    b = 1.0
    c = 1.0
    d = 1.0
    return c*exp(x) + b*exp(y) - d*x - a*y
end

function H_Q11f52_exp(x,y)
    a = Q11f52(1.0)
    b = Q11f52(1.0)
    c = Q11f52(1.0)
    d = Q11f52(1.0)
    return c*Q11f52(exp(x)) + b*Q11f52(exp(y)) - d*x - a*y
end

function H_Q4f59_exp(x,y)
    a = Q4f59(1.0)
    b = Q4f59(1.0)
    c = Q4f59(1.0)
    d = Q4f59(1.0)
    return c*Q4f59(exp(x)) + b*Q4f59(exp(y)) - d*x - a*y
end

ans_bigfloat = spm_lv_T([BigFloat("0.5"),BigFloat("0.25")],f_BigFloat_1,f_BigFloat_2,vec_spm_BigFloat,mat_spm_BigFloat,[BigFloat("1.0"),BigFloat("1.0"),BigFloat("1.0"),BigFloat("1.0")],BigFloat(2^(-7)),BigFloat,2^14,100,1.0e-20)
fig = Figure(size = (800,400))
x = exp.(ans_bigfloat[1])
y = exp.(ans_bigfloat[2])
ax1 = Axis(fig[1,1],xlabel = "t",title = "smp  BigFloat Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "x",ylabel = "y",title = "smp  BigFloat Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y) 
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_BigFloat_spm.pdf", fig)

fig = Figure(size = (600,400))
h = H_BigFloat_exp.(ans_bigfloat[1],ans_bigfloat[2]) .- H_BigFloat_exp(ans_bigfloat[1][1],ans_bigfloat[2][1])
ax = Axis(fig[1,1],xlabel = "t",title = "H error smp BigFloat Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_BigFloat_spm_h_error.pdf", fig)

ans_float64 = spm_lv_T([0.5,0.25],f_float64_1,f_float64_2,vec_spm_float64,mat_spm_float64,[1.0,1.0,1.0,1.0],2^(-7),Float64,2^(14),400,1.0e-20)
fig = Figure(size = (800,400))
x = exp.(ans_float64[1])
y = exp.(ans_float64[2])
ax1 = Axis(fig[1,1],xlabel = "t",title = "smp  Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "x",ylabel = "y")
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y) 
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_spm.pdf", fig)

fig = Figure(size = (600,400))
h = H_float64_exp.(ans_float64[1],ans_float64[2]) .- H_float64_exp(ans_float64[1][1],ans_float64[2][1])
ax = Axis(fig[1,1],xlabel = "t",title = "H error Float64 Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_spm_h_error.pdf", fig)

ans_float64_shift = spm_lv_T_shift([0.5,0.25],f_float64_1,f_float64_2,vec_spm_float64,mat_spm_float64,[1.0,1.0,1.0,1.0],2^(-7),[2.0,2.0],Float64,2^(14),400,1.0e-20)
fig = Figure(size = (800,400))
x = exp.(ans_float64_shift[1])
y = exp.(ans_float64_shift[2])
ax1 = Axis(fig[1,1],xlabel = "t",title = "smp Float64 shift Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "x",ylabel = "y")
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y) 
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_spm_shift.pdf", fig)

fig = Figure(size = (600,400))
h = H_float64_exp.(ans_float64_shift[1],ans_float64_shift[2]) .- H_float64_exp(ans_float64_shift[1][1],ans_float64_shift[2][1])
ax = Axis(fig[1,1],xlabel = "t",title = "H error Float64 shift Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_spm_shift_h_error.pdf", fig)

ans_Q11f52 = spm_lv_T([Q11f52(0.5),Q11f52(0.25)],f_Q11f52_1,f_Q11f52_2,vec_spm_Q11f52,mat_spm_Q11f52,[Q11f52(1.0),Q11f52(1.0),Q11f52(1.0),Q11f52(1.0)],Q11f52(2^(-7)),Q11f52,2^(14),400,1.0e-20)
fig = Figure(size = (800,400))
x = exp.(ans_Q11f52[1])
y = exp.(ans_Q11f52[2])
ax1 = Axis(fig[1,1],xlabel = "t",title = "smp test Q11f52")
ax2 = Axis(fig[1,2],xlabel = "x",ylabel = "y")
t = [i for i in 0:0.01:10]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y) 
fig

fig = Figure(size = (500,400))
h = H_Q11f52_exp.(ans_Q11f52[1],ans_Q11f52[2]) .- H_Q11f52_exp(ans_Q11f52[1][1],ans_Q11f52[2][1])
ax = Axis(fig[1,1],xlabel = "t",title = "H Q11f52 no shift")
t = [i for i in 0:0.01:10]
l = lines!(ax,t,h)
fig

ans_Q11f52_shift = spm_lv_T_shift([Q11f52(0.5),Q11f52(0.25)],f_Q11f52_1,f_Q11f52_2,vec_spm_Q11f52,mat_spm_Q11f52,[Q11f52(1.0),Q11f52(1.0),Q11f52(1.0),Q11f52(1.0)],Q11f52(2^-7),[Q11f52(3.0),Q11f52(3.0)],Q11f52,2^(14),400,1.0e-20)
fig = Figure(size = (800,400))
x = exp.(ans_Q11f52_shift[1])
y = exp.(ans_Q11f52_shift[2])
ax1 = Axis(fig[1,1],xlabel = "t",title = "smp  Q11f52 shift Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "x",ylabel = "y")
t = [i for i in 0:2^(-7):2^(7)]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y) 
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q11f52_spm_shift.pdf", fig)

fig = Figure(size = (600,400))
h = H_Q11f52_exp.(ans_Q11f52_shift[1],ans_Q11f52_shift[2]) .- H_Q11f52_exp(ans_Q11f52_shift[1][1],ans_Q11f52_shift[2][1])
ax = Axis(fig[1,1],xlabel = "t",title = "H error Q11f52 shift Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q11f52_spm_shift_h_error.pdf", fig)

ans_Q4f59_shift = spm_lv_T_shift([Q4f59(0.5),Q4f59(0.25)],f_Q4f59_1,f_Q4f59_2,vec_spm_Q4f59,mat_spm_Q4f59,[Q4f59(1.0),Q4f59(1.0),Q4f59(1.0),Q4f59(1.0)],Q4f59(2^-7),[Q4f59(1.0),Q4f59(1.0)],Q4f59,2^(14),400,1.0e-20)
fig = Figure(size = (800,400))
x = exp.(ans_Q4f59_shift[1])
y = exp.(ans_Q4f59_shift[2])
ax1 = Axis(fig[1,1],xlabel = "t",title = "smp Q4f59 shift Δt = 2^(-7)")
ax2 = Axis(fig[1,2],xlabel = "x",ylabel = "y")
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y) 
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q4f59_spm_shift.pdf", fig)

fig = Figure(size = (600,400))
h = H_Q4f59_exp.(ans_Q4f59_shift[1],ans_Q4f59_shift[2]) .- H_Q4f59_exp(ans_Q4f59_shift[1][1],ans_Q4f59_shift[2][1])
ax = Axis(fig[1,1],xlabel = "t",title = "H error Q4f59 shift Δt = 2^(-7)")
t = [i for i in 0:2^(-7):2^7]
l = lines!(ax,t,h)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_Q4f59_spm_shift_h_error.pdf", fig)

#それぞれの誤差を比較
fig = Figure(size = (600,400))
ax = Axis(fig[1,1],xlabel = "t",title = "H error")
h1 = H_float64_exp.(ans_float64[1],ans_float64[2]) .- H_float64_exp(ans_float64[1][1],ans_float64[2][1])
h2 = H_float64_exp.(ans_float64_shift[1],ans_float64_shift[2]) .- H_float64_exp(ans_float64_shift[1][1],ans_float64_shift[2][1])
h3 = H_Q11f52_exp.(ans_Q11f52_shift[1],ans_Q11f52_shift[2]) .- H_Q11f52_exp(ans_Q11f52_shift[1][1],ans_Q11f52_shift[2][1])
h4 = H_Q4f59_exp.(ans_Q4f59_shift[1],ans_Q4f59_shift[2]) .- H_Q4f59_exp(ans_Q4f59_shift[1][1],ans_Q4f59_shift[2][1])
t = [i for i in 0:2^(-7):2^7]
l1 = lines!(ax,t,h1,color = :red)
l2 = lines!(ax,t,h2,color = :blue)
l3 = lines!(ax,t,h3,color = cgrad(:Accent_4)[1])
l4 = lines!(ax,t,h4,color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,l2,l3,l4],["Float64","Float64 shift","Q11f52 shift","Q4f59 shift"],position = :rb)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_h_error.pdf", fig)

#解の精度の比較
fig = Figure(size = (600,400))
ax = Axis(fig[1,1],xlabel = "t",title = "smp Numerical error")
t = [i for i in 0:2^(-7):2^7]
err1 = err_norm_2d(ans_bigfloat[1],ans_bigfloat[2],ans_float64[1],ans_float64[2])
err2 = err_norm_2d(ans_bigfloat[1],ans_bigfloat[2],ans_float64_shift[1],ans_float64_shift[2])
err3 = err_norm_2d(ans_bigfloat[1],ans_bigfloat[2],ans_Q11f52_shift[1],ans_Q11f52_shift[2])
err4 = err_norm_2d(ans_bigfloat[1],ans_bigfloat[2],ans_Q4f59_shift[1],ans_Q4f59_shift[2])
l1 = lines!(ax,t,err1,color = :red)
l2 = lines!(ax,t,err2,color = :blue)
l3 = lines!(ax,t,err3,color = cgrad(:Accent_4)[1])
l4 = lines!(ax,t,err4,color = cgrad(:Accent_4)[2])
axislegend(ax,[l1,l2,l3,l4],["Float64","Float64 shift","Q11f52 shift","Q4f59 shift"],position = :rb)
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_smp_error.pdf", fig)