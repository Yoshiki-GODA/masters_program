using LinearAlgebra,FixedPointNumbers,CairoMakie
include("../numerical_scheme.jl")
include("../numerical_scheme3.jl")
setprecision(BigFloat,256)

#構造保存数値解法の実装
#Newton法
#=function my_newton(H,dH_inv,x_ini,stepsize,Type,itr,tol)
    T = Type
    x = x_ini
    
    for i in 1:itr
        if norm(H(x),2) < tol
            break
        end
        x = x - dH_inv(x)*H(x)
    end
    return x,norm(H(x))
end=#

#構造保存数値解法
#=function spm_lv_T(v_0,f_1,f_2,vec_t,mat_inv_t,p,stepsize,Type,itr,itr_newton,tol_newton)
    T = Type
    s = T(stepsize)
    result1 = T[log(v_0[1])]
    result2 = T[log(v_0[2])]
    sizehint!(result1,itr+1)
    sizehint!(result2,itr+1)
    S = T(0.0)
    for i in 1:itr
        v_ini = [euler_2(f_1,f_2,exp(result1[end]),exp(result2[end]),s,T,1)[1][end],euler_2(f_1,f_2,result1[end],result2[end],s,T,1)[2][end]]
        vec(v_new) = vec_t(v_new,[result1[end],result2[end]],s,p)
        mat(v_new) = mat_inv_t(v_new,[result1[end],result2[end]],s,p)
        v = my_newton(vec,mat,v_ini,s,T,itr_newton,tol_newton)[1]
        push!(result1,v[1])
        push!(result2,v[2])
    end
    return result1,result2
end=#

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

#更新式の行列
function vec_spm_float64(v_new::Vector{Float64},v_old::Vector{Float64},s::Float64,p::Vector{Float64})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ]
end

function vec_spm_BigFloat(v_new::Vector{BigFloat},v_old::Vector{BigFloat},s::BigFloat,p::Array{BigFloat,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ]
end

function vec_spm_Q11f52(v_new::Vector{Q11f52},v_old::Vector{Q11f52},s::Q11f52,p::Array{Q11f52,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ]
end

function vec_spm_Q3f60(v_new::Vector{Q3f60},v_old::Vector{Q3f60},s::Q3f60,p::Array{Q3f60,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ] 
end

#更新式の行列のヤコビアンの逆行列
function mat_inv_spm_float64(v_new::Vector{Float64},v_old::Vector{Float64},s::Float64,p::Vector{Float64})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    a1 = s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2)))
    a2 = s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2)))
    b1 = 1 - a1*a2
    c1 = 1/b1
    c2 = -a1/b1
    c3 = -a2/b1
    return [
        c1 c2
        c3 c1
    ]
end

function mat_inv_spm_BigFloat(v_new::Vector{BigFloat},v_old::Vector{BigFloat},s::BigFloat,p::Array{BigFloat,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    a1 = s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2)))
    a2 = s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2)))
    b1 = 1 - a1*a2
    c1 = 1/b1
    c2 = -a1/b1
    c3 = -a2/b1
    return [
        c1 c2
        c3 c1
    ]
end

function mat_inv_spm_Q11f52(v_new::Array{Q11f52,1},v_old::Array{Q11f52,1},s::Q11f52,p::Array{Q11f52,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    a1 = Q11f52(s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2))))
    a2 = Q11f52(s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2))))
    b1 = 1 - a1*a2
    c1 = 1/b1
    c2 = -a1/b1
    c3 = -a2/b1
    return [
        c1 c2
        c3 c1
    ]
end

function mat_inv_spm_Q3f60(v_new::Array{Q3f60,1},v_old::Array{Q3f60,1},s::Q3f60,p::Array{Q3f60,1})
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new 
    a1 = Q3f60(s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2))))
    a2 = Q3f60(s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2))))
    b1 = 1 - a1*a2
    c1 = 1/b1
    c2 = -a1/b1
    c3 = -a2/b1
    return [
        c1 c2
        c3 c1
    ]
end

#保存量の計算
function H_float64_exp(x,y)
    a = 1.0
    b = 1.0
    c = 1.0
    d = 1.0
    return c*exp(x) + b*exp(y) - d*x - a*y
end

function H_BigFloat_exp(x,y)
    a = BigFloat("1.0")
    b = BigFloat("1.0")
    c = BigFloat("1.0")
    d = BigFloat("1.0")
    return c*exp(x) + b*exp(y) - d*x - a*y
end

function H_Q11f52_exp(x,y)
    a = Q11f52(1.0)
    b = Q11f52(1.0)
    c = Q11f52(1.0)
    d = Q11f52(1.0)
    return c*Q11f52(exp(x)) + b*Q11f52(exp(y)) - d*x - a*y
end

function H_Q3f60_exp(x,y)
    a = Q3f60(1.0)
    b = Q3f60(1.0)
    c = Q3f60(1.0)
    d = Q3f60(1.0)
    return c*Q3f60(exp(x)) + b*Q3f60(exp(y)) - d*x - a*y
end

#spm
fig = Figure(size =(800, 400))
ax1 = Axis(fig[1,1], xlabel = "t",title = "Lotka-Volterra SPM Float64 Δt = 2^(-7)")
ax2 = Axis(fig[1,2], xlabel = "x",ylabel = "y",title = "Lotka-Volterra SPM Float64 Δt = 2^(-7)")
x = spm_lv_T([0.5,0.25],f_float64_1,f_float64_2,vec_spm_float64,mat_inv_spm_float64,[1.0,1,0,1.0,1.0],2^(-7),Float64,Int(2^(14)),1000,1e-20)
y = spm_lv_T([0.5,0.25],f_float64_1,f_float64_2,vec_spm_float64,mat_inv_spm_float64,[1.0,1,0,1.0,1.0],2^(-7),Float64,Int(2^(9)),1000,1e-20)[2]
t = [i for i in 0:2^(-7):2^2]
l1 = lines!(ax1,t,x)
l2 = lines!(ax1,t,y)
l3 = lines!(ax2,x,y)
axislegend(ax1,[l1,l2],["x","y"])
fig
save("lotka-volterra/fig_lotkavolterra/lotkavolterra_2^(-7)_float64_spm.pdf", fig)

fig = Figure(size =(800, 400))