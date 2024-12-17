using Plots
#Lorenz方程式の定数
σ_0 = 10.0
ρ_0 = 28.0
β_0 = 2.6640625

#刻み幅
h_0 = 0.001

using FixedPointNumbers

#初期値
x_0 = 50.0
y_0 = 50.0
z_0 = 50.0

#オイラー法
function  euler(x::T,y::T,z::T,σ::T,ρ::T,β::T,h::T) where T
    x_new = x + h*(σ*(y-x))
    y_new = y + h*(x*(ρ-z)-y)
    z_new = z + h*(x*y - β*z)
    return x_new, y_new, z_new
end

#BigFloatで実行
σ_big = BigFloat(σ_0)
ρ_big = BigFloat(ρ_0)
β_big = BigFloat(β_0) 

h_big = BigFloat(h_0)

x_big = BigFloat(x_0)
y_big = BigFloat(y_0)
z_big = BigFloat(z_0)

ans2 = euler(x_big,y_big,z_big,σ_big,ρ_big,β_big,h_big)

typeof(ans2)

x,y,z = x_big, y_big, z_big
sq_big = [x_big y_big z_big]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_big,ρ_big,β_big,h_big)
    sq_big = vcat(sq_big, [x y z])
end

sq2_big = [x_big y_big z_big]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_big,ρ_big,β_big,h_big)
    sq2_big = vcat(sq2_big, [x y z])
end

typeof(sq_big[:,1])

#計算結果を図示
plot(sq2_big[:,1],sq2_big[:,2],sq2_big[:,3],title = "lorenz bigfloat")

maximum(sq_big[:,1])
minimum(sq_big[:,1])

maximum(sq_big[:,2])
minimum(sq_big[:,2])

maximum(sq_big[:,3])
minimum(sq_big[:,3])

#単精度浮動小数点数で実行
σ_float32 = Float32(σ_0)
ρ_float32 = Float32(ρ_0)
β_float32 = Float32(β_0)

h_float32 = Float32(h_0)

x_float32 = Float32(x_0)
y_float32 = Float32(y_0)
z_float32 = Float32(z_0)

ans1 = euler(x_float32,y_float32,z_float32,σ_float32,ρ_float32,β_float32,h_float32)

typeof(ans1)

x,y,z = x_float32, y_float32, z_float32
sq_float32 = [x_float32 y_float32 z_float32]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_float32,ρ_float32,β_float32,h_float32)
    sq_float32 = vcat(sq_float32, [x y z])
end

#計算結果を図示
plot(sq_float32[:,1],sq_float32[:,2],sq_float32[:,3],title = "lorenz float32")

err_big_float32_sq = sq_big .- sq_float32

plot(err_big_float32_sq[:,1],title = "x_error bigfloat vs float32")
plot(err_big_float32_sq[:,2],title = "y_error bigfloat vs float32")
plot(err_big_float32_sq[:,3],title = "z_error bigfloat vs float32")

#倍精度浮動小数点で実行
σ_float64 = Float64(σ_0)
ρ_float64 = Float64(ρ_0)
β_float64 = Float64(β_0)

h_float64 = Float64(h_0)

x_float64 = Float64(x_0)
y_float64 = Float64(y_0)
z_float64 = Float64(z_0)

ans1 = euler(x_float64,y_float64,z_float64,σ_float64,ρ_float64,β_float64,h_float64)

typeof(ans1)

x,y,z = x_float64, y_float64, z_float64
sq_float64 = [x_float64 y_float64 z_float64]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_float64,ρ_float64,β_float64,h_float64)
    sq_float64 = vcat(sq_float64, [x y z])
end

typeof(sq_float64)

#計算結果を図示
plot(sq_float64[:,1],sq_float64[:,2],sq_float64[:,3],title = "lorenz float64")

#bigfoaltとの比較
err_big_float64_sq = sq_big .- sq_float64

iter_sq = [j for j in 1:1000001]

maximum(err_big_float64_sq[:,1])

plot(err_big_float64_sq[:,1],title = "x_error bigfloat vs float64")
plot(err_big_float64_sq[:,2],title = "y_error bigfloat vs float64")
plot(err_big_float64_sq[:,3],title = "z_error bigfloat vs float64")

#64ビットの固定小数点数で実行
Fixed{Int64,32}

σ_fixed64 = Q31f32(σ_0)
ρ_fixed64 = Q31f32(ρ_0)
β_fixed64 = Q31f32(β_0)

h_fixed64 = Q31f32(h_0) 

x_fixed64 = Q31f32(x_0)
y_fixed64 = Q31f32(y_0)
z_fixed64 = Q31f32(z_0)

ans3 = euler(x_fixed64,y_fixed64,z_fixed64,σ_fixed64,ρ_fixed64,β_fixed64,h_fixed64)

typeof(ans3)

x,y,z = x_fixed64, y_fixed64, z_fixed64
sq_fixed64 = [x_fixed64 y_fixed64 z_fixed64]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_fixed64,ρ_fixed64,β_fixed64,h_fixed64)
    sq_fixed64 = vcat(sq_fixed64, [x y z])
end

#計算結果を図示
plot(sq_fixed64[:,1],sq_fixed64[:,2],sq_fixed64[:,3],title = "lorenz Q31f32")

err_big_Q31f32_sq = sq_big .- sq_fixed64

plot(err_big_Q31f32_sq[:,1],title = "x_error bigfloat vs Q31f32")
plot(err_big_Q31f32_sq[:,2],title = "y_error bigfloat vs Q31f32")
plot(err_big_Q31f32_sq[:,3],title = "z_error bigfloat vs Q31f32")

maximum(sq_fixed64[:,1])

#64ビットの固定小数点数で実行
Fixed{Int64,42}

σ_Q21f42 = Q21f42(σ_0)
ρ_Q21f42 = Q21f42(ρ_0)
β_Q21f42 = Q21f42(β_0)

h_Q21f42 = Q21f42(h_0)

x_Q21f42 = Q21f42(x_0)
y_Q21f42 = Q21f42(y_0)
z_Q21f42 = Q21f42(z_0)  

ans4 = euler(x_Q21f42,y_Q21f42,z_Q21f42,σ_Q21f42,ρ_Q21f42,β_Q21f42,h_Q21f42)

typeof(ans4)

x,y,z = x_Q21f42, y_Q21f42, z_Q21f42
sq_Q21f42 = [x_Q21f42 y_Q21f42 z_Q21f42]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_Q21f42,ρ_Q21f42,β_Q21f42,h_Q21f42)
    sq_Q21f42 = vcat(sq_Q21f42, [x y z])
end

#計算結果を図示
plot(sq_Q21f42[:,1],sq_Q21f42[:,2],sq_Q21f42[:,3],title = "lorenz Q21f42")

err_big_Q21f42_sq = sq_big .- sq_Q21f42

plot(err_big_Q21f42_sq[:,1],title = "x_error bigfloat vs Q21f42")
plot(err_big_Q21f42_sq[:,2],title = "y_error bigfloat vs Q21f42")
plot(err_big_Q21f42_sq[:,3],title = "z_error bigfloat vs Q21f42")

#64ビットの固定小数点数で実行
Fixed{Int64,52}

σ_Q11f52 = Q11f52(σ_0)
ρ_Q11f52 = Q11f52(ρ_0)
β_Q11f52 = Q11f52(β_0)

h_Q11f52 = Q11f52(h_0)

x_Q11f52 = Q11f52(x_0)
y_Q11f52 = Q11f52(y_0)
z_Q11f52 = Q11f52(z_0)

ans5 = euler(x_Q11f52,y_Q11f52,z_Q11f52,σ_Q11f52,ρ_Q11f52,β_Q11f52,h_Q11f52)

typeof(ans5)

x,y,z = x_Q11f52, y_Q11f52, z_Q11f52
sq_Q11f52 = [x_Q11f52 y_Q11f52 z_Q11f52]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_Q11f52,ρ_Q11f52,β_Q11f52,h_Q11f52)
    sq_Q11f52 = vcat(sq_Q11f52, [x y z])
end

#計算結果を図示
plot(sq_Q11f52[:,1],sq_Q11f52[:,2],sq_Q11f52[:,3],title = "lorenz Q11f52")

err_big_Q11f52_sq = sq_big .- sq_Q11f52

plot(err_big_Q11f52_sq[:,1],title = "x_error bigfloat vs Q11f52")
plot(err_big_Q11f52_sq[:,2],title = "y_error bigfloat vs Q11f52")
plot(err_big_Q11f52_sq[:,3],title = "z_error bigfloat vs Q11f52")

#64ビットの固定小数点数で実行
Fixed{Int64,56}

σ_Q7f56 = Q7f56(σ_0)
ρ_Q7f56 = Q7f56(ρ_0)
β_Q7f56 = Q7f56(β_0)

h_Q7f56 = Q7f56(h_0)

x_Q7f56 = Q7f56(x_0)
y_Q7f56 = Q7f56(y_0)
z_Q7f56 = Q7f56(z_0)

ans6 = euler(x_Q7f56,y_Q7f56,z_Q7f56,σ_Q7f56,ρ_Q7f56,β_Q7f56,h_Q7f56)

typeof(ans6)

x,y,z = x_Q7f56, y_Q7f56, z_Q7f56
sq_Q7f56 = [x_Q7f56 y_Q7f56 z_Q7f56]
for j in 1:100000
    x, y, z = euler(x,y,z,σ_Q7f56,ρ_Q7f56,β_Q7f56,h_Q7f56)
    sq_Q7f56 = vcat(sq_Q7f56, [x y z])
end

#計算結果を図示
plot(sq_Q7f56[:,1],sq_Q7f56[:,2],sq_Q7f56[:,3],title = "lorenz Q7f56") #失敗していた．．．

t = [i for i in 1:1e-4:100.0]