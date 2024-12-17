using GLMakie
using FixedPointNumbers
using DifferentialEquations
include("../numerical_scheme.jl")

#Lorenz方程式

#倍精度浮動小数点数で実行
#=σ = 10.0
ρ = 28.0
β = 8/3=#

#刻み幅1e-3

t_sq2 = range(0.0,100.0,step = 0.001)

function f1_fl64(x::Float64,y::Float64,z::Float64,t::Float64)
    σ = 10.0
    return σ*(y - x)
end

function f2_fl64(x::Float64,y::Float64,z::Float64,t::Float64)
    ρ = 28.0
    return x*(ρ - z) - y
end

function f3_fl64(x::Float64,y::Float64,z::Float64,t::Float64)
    β = 8/3
    return x*y - β*z
end


lorenz_float64_euler = euler_3(f1_fl64,f2_fl64,f3_fl64,0.1,0.1,0.1,0.001,Float64,100000)

plot(lorenz_float64_euler[1],lorenz_float64_euler[2],lorenz_float64_euler[3],title = "lorenz Float64 Euler",xlabel = "x",ylabel = "y",zlabel = "z",label = false)
savefig("lorenz_Float64_Euler.pdf")

#64ビットの固定小数点数で実行

#32bitを小数部分にする
Fixed{Int64,32}
Fixed{Int64,22}

σ_Q31f32 = Q31f32(σ)
ρ_Q31f32 = Q31f32(ρ)
β_Q31f32 = Q31f32(β)

function f_31f32_1(x::Q31f32,y::Q31f32,z::Q31f32,t::Q31f32)
    return σ_Q31f32*(y - x)
end

function f_31f32_2(x::Q31f32,y::Q31f32,z::Q31f32,t::Q31f32)
    return x*(ρ_Q31f32 - z) - y
end
function f_31f32_3(x::Q31f32,y::Q31f32,z::Q31f32,t::Q31f32)
    return x*y - β_Q31f32*z
end

#42bitを小数部分にする
σ_Q21f42 = Q21f42(σ)
ρ_Q21f42 = Q21f42(ρ)
β_Q21f42 = Q21f42(β)

function f_Q21f42_1(x::Q21f42,y::Q21f42,z::Q21f42,t::Q21f42)
    return σ_Q21f42*(y - x)
end
function f_Q21f42_2(x::Q21f42,y::Q21f42,z::Q21f42,t::Q21f42)
    return x*(ρ_Q21f42 - z) - y
end
function f_Q21f42_3(x::Q21f42,y::Q21f42,z::Q21f42,t::Q21f42)
    return x*y - β_Q21f42*z
end

lorenz_Q21f42_euler = euler_3(f_Q21f42_1,f_Q21f42_2,f_Q21f42_3,Q21f42(0.1),Q21f42(0.1),Q21f42(0.1),Q21f42(0.001),Q21f42,100000)
plot(lorenz_Q21f42_euler[1],lorenz_Q21f42_euler[2],lorenz_Q21f42_euler[3],title = "lorenz Q21f42 Euler",xlabel = "x",ylabel = "y",zlabel = "z",label = false)
savefig("lorenz_Q21f42_Euler.pdf")

lorenz_31f32_euler = euler_3(f_31f32_1,f_31f32_2,f_31f32_3,Q31f32(0.1),Q31f32(0.1),Q31f32(0.1),Q31f32(0.001),Q31f32,100000)
plot(lorenz_31f32_euler[1],lorenz_31f32_euler[2],lorenz_31f32_euler[3],title = "lorenz Q31f32 Euler")
savefig("lorenz_Q31f32_Euler.pdf")

lorenz_31f32_euler[1]

#bigfloatで実行
σ_big = big(σ)
ρ_big = big(ρ)
β_big = big(β)

function f_big_1(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    return σ_big*(y - x)
end
function f_big_2(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    return x*(ρ_big - z) - y
end
function f_big_3(x::BigFloat,y::BigFloat,z::BigFloat,t::BigFloat)
    return x*y - β_big*z
end

lorenz_big_euler = euler_3(f_big_1,f_big_2,f_big_3,big(0.1),big(0.1),big(0.1),big(0.001),BigFloat,100000)

lorenz_big_euler[1]

plot(lorenz_big_euler[1],lorenz_big_euler[2],lorenz_big_euler[3],title = "lorenz bigfloat Euler")
savefig("lorenz_bigfloat_Euler.pdf")

#各成分ごとの誤差を計算
#float64とbigfloat
err_float64_vs_big_euler_1 = lorenz_float64_euler[1] .- lorenz_big_euler[1]
plot(err_float64_vs_big_euler_1,title = "x_error float64 vs bigfloat")

err_float64_vs_big_euler_2 = lorenz_float64_euler[2] .- lorenz_big_euler[2]
plot(err_float64_vs_big_euler_2,title = "y_error float64 vs bigfloat")

err_float64_vs_big_euler_3 = lorenz_float64_euler[3] .- lorenz_big_euler[3]
plot(err_float64_vs_big_euler_3,title = "z_error float64 vs bigfloat")

#Q31f32とbigfloat
err_Q31f32_vs_big_euler_1 = lorenz_31f32_euler[1] .- lorenz_big_euler[1]
plot(err_Q31f32_vs_big_euler_1,title = "x_error Q31f32 vs bigfloat")

err_Q31f32_vs_big_euler_2 = lorenz_31f32_euler[2] .- lorenz_big_euler[2]
plot(err_Q31f32_vs_big_euler_2,title = "y_error Q31f32 vs bigfloat")

err_Q31f32_vs_big_euler_3 = lorenz_31f32_euler[3] .- lorenz_big_euler[3]
plot(err_Q31f32_vs_big_euler_3,title = "z_error Q31f32 vs bigfloat")

lorenz_31f32_euler[1]
lorenz_big_euler[1]

#float64とQ21f42
err_Q21f42_vs_big_1 = lorenz_Q21f42_euler[1] .- lorenz_big_euler[1]
plot(err_Q21f42_vs_big_1,title = "x_error Q21f42 vs bigfloat")

err_Q21f42_vs_big_2 = lorenz_Q21f42_euler[2] .- lorenz_big_euler[2]
plot(err_Q21f42_vs_big_2,title = "y_error Q21f42 vs bigfloat")

err_Q21f42_vs_big_3 = lorenz_Q21f42_euler[3] .- lorenz_big_euler[3]
plot(err_Q21f42_vs_big_3,title = "z_error Q21f42 vs bigfloat")

#ノルムでの誤差を計算
#float64とbigfloat
err__float64_vs_big_euler_norm = []
for i in 1:length(lorenz_float64_euler[1])
    x = err_float64_vs_big_euler_1[i]
    y = err_float64_vs_big_euler_2[i]
    z = err_float64_vs_big_euler_3[i]
    push!(err__float64_vs_big_euler_norm, sqrt(x^2+y^2+z^2))
end

plot(err__float64_vs_big_euler_norm,title = "norm_error_float64_vs_bigfloat_euler")
savefig("norm_error_float64_vs_bigfloat_euler.pdf")

#Q31f32とbigfloat
err_Q31f32_vs_big_euler_norm = []
for i in 1:length(lorenz_31f32_euler[1])
    x = err_Q31f32_vs_big_euler_1[i]
    y = err_Q31f32_vs_big_euler_2[i]
    z = err_Q31f32_vs_big_euler_3[i]
    push!(err_Q31f32_vs_big_euler_norm, sqrt(x^2+y^2+z^2))
end

plot(err_Q31f32_vs_big_euler_norm,title = "norm_error_Q31f32_vs_bigfloat_euler")
savefig("norm_error_Q31f32_vs_bigfloat_euler.pdf")

#Q21f42とbigfloat
err_Q21f42_vs_big_euler_norm = []
for i in 1:length(lorenz_Q21f42_euler[1])
    x = err_Q21f42_vs_big_1[i]
    y = err_Q21f42_vs_big_2[i]
    z = err_Q21f42_vs_big_3[i]
    push!(err_Q21f42_vs_big_euler_norm, sqrt(x^2+y^2+z^2))
end

plot(err_Q21f42_vs_big_euler_norm,title = "norm_error_Q21f42_vs_bigfloat_euler")
savefig("norm_error_Q21f42_vs_bigfloat_euler.pdf")

#誤差を比較
plot(err__float64_vs_big_euler_norm,label = "float64")
plot!(err_Q31f32_vs_big_euler_norm,label = "Q31f32")
plot!(err_Q21f42_vs_big_euler_norm,label = "Q21f42")
savefig("error_compare_lorenz_euler.pdf")

#誤差の最大値を比較
maximum(err__float64_vs_big_euler_norm)
maximum(err_Q31f32_vs_big_euler_norm)
maximum(err_Q21f42_vs_big_euler_norm)


#ルンゲ・クッタ法
#float64
lorenz_float64_rk4 = rk4_3(f1_fl64,f2_fl64,f3_fl64,0.1,0.1,0.1,0.001,Float64,100000)

plot(lorenz_float64_rk4[1],lorenz_float64_rk4[2],lorenz_float64_rk4[3],title = "lorenz float64 rk4")
savefig("lorenz_float64_rk4.pdf")

#Q31f32
lorenz_Q31f32_rk4 = rk4_3(f_31f32_1,f_31f32_2,f_31f32_3,Q31f32(0.1),Q31f32(0.1),Q31f32(0.1),Q31f32(0.001),Q31f32,100000)

plot(lorenz_Q31f32_rk4[1],lorenz_Q31f32_rk4[2],lorenz_Q31f32_rk4[3],title = "lorenz Q31f32 rk4")
savefig("lorenz_Q31f32_rk4.pdf")

#Q21f42
lorenz_Q21f42_rk4 = rk4_3(f_Q21f42_1,f_Q21f42_2,f_Q21f42_3,Q21f42(0.1),Q21f42(0.1),Q21f42(0.1),Q21f42(0.001),Q21f42,100000)

plot(lorenz_Q21f42_rk4[1],lorenz_Q21f42_rk4[2],lorenz_Q21f42_rk4[3],title = "lorenz Q21f42 rk4")
savefig("lorenz_Q21f42_rk4.pdf")

#bigfloatで実行
lorenz_big_rk4 = rk4_3(f_big_1,f_big_2,f_big_3,big(0.1),big(0.1),big(0.1),big(0.001),BigFloat,100000)

plot(lorenz_big_rk4[1],lorenz_big_rk4[2],lorenz_big_rk4[3],title = "lorenz bigfloat rk4")
savefig("lorenz_big_rk4.pdf")

#各成分ごとの誤差を計算
#float64とbigfloat
err_float64_vs_big_rk4_1 = lorenz_big_rk4[1] .- lorenz_float64_rk4[1]
err_float64_vs_big_rk4_2 = lorenz_big_rk4[2] .- lorenz_float64_rk4[2]
err_float64_vs_big_rk4_3 = lorenz_big_rk4[3] .- lorenz_float64_rk4[3]

#Q31f32とbigfloat
err_Q31f32_vs_big_rk4_1 = lorenz_big_rk4[1] .- lorenz_Q31f32_rk4[1]
err_Q31f32_vs_big_rk4_2 = lorenz_big_rk4[2] .- lorenz_Q31f32_rk4[2]
err_Q31f32_vs_big_rk4_3 = lorenz_big_rk4[3] .- lorenz_Q31f32_rk4[3]

#Q21f42とbigfloat
err_Q21f42_vs_big_rk4_1 = lorenz_big_rk4[1] .- lorenz_Q21f42_rk4[1]
err_Q21f42_vs_big_rk4_2 = lorenz_big_rk4[2] .- lorenz_Q21f42_rk4[2]
err_Q21f42_vs_big_rk4_3 = lorenz_big_rk4[3] .- lorenz_Q21f42_rk4[3]

#ノルムでの誤差を計算
#float64とbigfloat
err_float64_vs_big_rk_norm = []
for i in 1:length(lorenz_float64_rk4[1])
    x = err_float64_vs_big_rk4_1[i]
    y = err_float64_vs_big_rk4_2[i]
    z = err_float64_vs_big_rk4_3[i]
    push!(err_float64_vs_big_rk_norm, sqrt(x^2+y^2+z^2))
end

plot(err_float64_vs_big_rk_norm,title = "norm_error float64 vs bigfloat rk")
savefig("norm_error_float64_vs_bigfloat_rk.pdf")

#Q31f32とbigfloat
err_Q31f32_vs_big_rk_norm = []
for i in 1:length(lorenz_Q31f32_rk4[1])
    x = err_Q31f32_vs_big_rk4_1[i]
    y = err_Q31f32_vs_big_rk4_2[i]
    z = err_Q31f32_vs_big_rk4_3[i]
    push!(err_Q31f32_vs_big_rk_norm, sqrt(x^2+y^2+z^2))
end

plot(err_Q31f32_vs_big_rk_norm,title = "norm_error Q31f32 vs bigfloat rk")
savefig("norm_error_Q31f32_vs_bigfloat_rk.pdf")

#Q21f42とbigfloat
err_Q21f42_vs_big_rk_norm = []
for i in 1:length(lorenz_Q21f42_rk4[1])
    x = err_Q21f42_vs_big_rk4_1[i]
    y = err_Q21f42_vs_big_rk4_2[i]
    z = err_Q21f42_vs_big_rk4_3[i]
    push!(err_Q21f42_vs_big_rk_norm, sqrt(x^2+y^2+z^2))
end

plot(err_Q21f42_vs_big_rk_norm,title = "norm_error Q21f42 vs bigfloat rk")
savefig("norm_error_Q21f42_vs_bigfloat_rk.pdf")

#誤差を比較
plot(err_float64_vs_big_rk_norm,label = "norm_error float64 vs bigfloat rk")
plot!(err_Q31f32_vs_big_rk_norm,label = "norm_error Q31f32 vs bigfloat rk")
plot!(err_Q21f42_vs_big_rk_norm,label = "norm_error Q21f42 vs bigfloat rk")  
savefig("norm_error_compare_rk.pdf")

#誤差の最大値を比較
maximum(err_float64_vs_big_rk_norm)
maximum(err_Q31f32_vs_big_rk_norm)
maximum(err_Q21f42_vs_big_rk_norm)

#微分方程式ソルバーを実行
function Lorenz!(du,u,p,t)
    σ,ρ,β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

tspan = (big(0.0),big(100.0))
u0 = [big(0.1);big(0.1);big(0.1)]
p = [big(10.0),big(28.0),big(8/3)]

prob = ODEProblem(Lorenz!,u0,tspan,p)
sol = solve(prob,Vern9(),reltol=1e-8,abstol=1e-8,saveat = big(0.001))
typeof(sol[1,:])
plot(sol[1,:],sol[2,:],sol[3,:],xlabel = "x",ylabel = "y",zlabel = "z")

#bigfloatでの計算結果を比較
err_euler_true_sq1 = sol[1,:] .- lorenz_big_euler[1]
err_euler_true_sq2 = sol[2,:] .- lorenz_big_euler[2]
err_euler_true_sq3 = sol[3,:] .- lorenz_big_euler[3]

#ノルムでの誤差を計算
err_euler_true_norm = []
for i in 1:length(lorenz_big_euler[1])
    x = err_euler_true_sq1[i]
    y = err_euler_true_sq2[i]
    z = err_euler_true_sq3[i]
    push!(err_euler_true_norm, sqrt(x^2+y^2+z^2))
end

plot(t_sq,err_euler_true_norm,yscale =:log10,xlabel = "t",ylabel = "norm_error",ylims = (1e-8,1e3),ytick = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3],title = "norm_error bigfloat_euler",label = "norm_error")

err_rk4_ture_sq1 = sol[1,:] .- lorenz_big_rk4[1]
err_rk4_ture_sq2 = sol[2,:] .- lorenz_big_rk4[2]
err_rk4_ture_sq3 = sol[3,:] .- lorenz_big_rk4[3]

#ノルムでの誤差を計算
err_rk4_ture_norm = []
for i in 1:length(lorenz_big_rk4[1])
    x = err_rk4_ture_sq1[i] 
    y = err_rk4_ture_sq2[i]
    z = err_rk4_ture_sq3[i]
    push!(err_rk4_ture_norm, sqrt(x^2+y^2+z^2))
end

plot(t_sq,err_rk4_ture_norm,yscale =:log10,xlabel = "t",ylabel = "norm_error",ylims = (1e-8,1e3),ytick = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3],title = "norm_error bigfloat_rk4",label = "rk4")
plot!(t_sq,err_euler_true_norm,label = "euler")

#float64での数値解と厳密解を比較
err_euler_fl64_ture_sq1 = sol[1,:] .- lorenz_float64_euler[1]
err_euler_fl64_ture_sq2 = sol[2,:] .- lorenz_float64_euler[2]
err_euler_fl64_ture_sq3 = sol[3,:] .- lorenz_float64_euler[3]

err_rk4_fl64_ture_sq1 = sol[1,:] .- lorenz_float64_rk4[1]
err_rk4_fl64_ture_sq2 = sol[2,:] .- lorenz_float64_rk4[2]   
err_rk4_fl64_ture_sq3 = sol[3,:] .- lorenz_float64_rk4[3]

#ノルムでの誤差を計算
err_euler_fl64_ture_norm = []
for i in 1:length(lorenz_float64_euler[1])
    x = err_euler_fl64_ture_sq1[i]  
    y = err_euler_fl64_ture_sq2[i]
    z = err_euler_fl64_ture_sq3[i]
    push!(err_euler_fl64_ture_norm, sqrt(x^2+y^2+z^2))
end 

err_rk4_fl64_ture_norm = []
for i in 1:length(lorenz_float64_rk4[1])
    x = err_rk4_fl64_ture_sq1[i] 
    y = err_rk4_fl64_ture_sq2[i]
    z = err_rk4_fl64_ture_sq3[i]
    push!(err_rk4_fl64_ture_norm, sqrt(x^2+y^2+z^2))
end

plot(t_sq,err_euler_fl64_ture_norm,yscale = :log10, xlabel = "t",ylabel = "norm_error",ylims = (1e-8,1e3),ytick = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3],title = "norm_error float64_euler",label = "float64 euler")
plot!(t_sq,err_rk4_ture_norm,label = "rk4")

#Q21f42での数値解と厳密解を比較
err_euler_Q21f42_ture_sq1 = sol[1,:] .- lorenz_Q21f42_euler[1]
err_euler_Q21f42_ture_sq2 = sol[2,:] .- lorenz_Q21f42_euler[2]
err_euler_Q21f42_ture_sq3 = sol[3,:] .- lorenz_Q21f42_euler[3]

err_rk4_Q21f42_ture_sq1 = sol[1,:] .- lorenz_Q21f42_rk4[1]
err_rk4_Q21f42_ture_sq2 = sol[2,:] .- lorenz_Q21f42_rk4[2]   
err_rk4_Q21f42_ture_sq3 = sol[3,:] .- lorenz_Q21f42_rk4[3]

#ノルムでの誤差を計算
err_euler_Q21f42_ture_norm = []
for i in 1:length(lorenz_Q21f42_euler[1])
    x = err_euler_Q21f42_ture_sq1[i]  
    y = err_euler_Q21f42_ture_sq2[i]
    z = err_euler_Q21f42_ture_sq3[i]
    push!(err_euler_Q21f42_ture_norm, sqrt(x^2+y^2+z^2))
end

err_rk4_Q21f42_ture_norm = []
for i in 1:length(lorenz_Q21f42_rk4[1])
    x = err_rk4_Q21f42_ture_sq1[i] 
    y = err_rk4_Q21f42_ture_sq2[i]
    z = err_rk4_Q21f42_ture_sq3[i]
    push!(err_rk4_Q21f42_ture_norm, sqrt(x^2+y^2+z^2))
end

plot(t_sq,err_euler_Q21f42_ture_norm,yscale = :log10, xlabel = "t",ylabel = "norm_error",ylims = (1e-8,1e3),ytick = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3],title = "norm_error Q21f42",label = "Q21f42 euler")
plot!(t_sq,err_rk4_Q21f42_ture_norm,label = "rk4")

#同じ解法ごとの数値解と厳密解を比較
plot(t_sq,err_euler_true_norm,yscale = :log10, xlabel = "t",ylabel = "norm_error",ylims = (1e-8,1e3),ytick = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3],title = "norm_error_eulur",label = "ture_err")
plot!(t_sq,err_euler_fl64_ture_norm,label = "float64")
plot!(t_sq,err_euler_Q21f42_ture_norm,label = "Q21f42")

plot(t_sq,err_rk4_ture_norm,yscale = :log10, xlabel = "t",ylabel = "norm_error",ylims = (1e-8,1e3),ytick = [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3],title = "norm_error_rk",label = "true_err")
plot!(t_sq,err_rk4_fl64_ture_norm,label = "float64")
plot!(t_sq,err_rk4_Q21f42_ture_norm,label = "Q21f42")