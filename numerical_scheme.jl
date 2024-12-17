using FixedPointNumbers

#誤差の計算
function err(x_ex,x)
    x_err = x_ex .- x
end

function err_norm_2d(x_ex,y_ex,x,y)
    x_err = x_ex .- x
    y_err = y_ex .- y
    return sqrt.((x_err).^2 .+ (y_err).^2)
end

function err_norm_3d(x_ex,y_ex,z_ex,x,y,z)
    x_err = x_ex .- x
    y_err = y_ex .- y
    z_err = z_ex .- z
    return sqrt.((x_err).^2 .+ (y_err).^2 .+ (z_err).^2)
end

#生成データの削減
function reduce_data(data,cutstep)
    n = cutstep
    T = typeof(data[1])
    cut_data = T[]
    for i in 1:n:length(data)
        push!(cut_data,data[i])
    end
    return cut_data
end

#オイラー法
function euler_1(f1,x0,stepsize,T,N) 
    s = T(stepsize)
    result1 = T[x0]
    sizehint!(result1,N+1)
    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        S+=s
        x_new = x1 + s*f1(x1,S)
        push!(result1,x_new)
    end
    return result1
end

function euler_2(f1,f2,x0,y0,stepsize,T,N)
    s = T(stepsize)
    result1 = T[x0]
    result2 = T[y0]
    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]
        S+=s
        x_new = x1 + s*f1(x1,y1,S)
        y_new = y1 + s*f2(x1,y1,S)
        push!(result1,x_new)
        push!(result2,y_new)
    end
    return result1, result2
end

function euler_3(f1,f2,f3,x0,y0,z0,stepsize,T,N)
    s = T(stepsize)
    result1 = T[x0]
    result2 = T[y0]
    result3 = T[z0]
    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    sizehint!(result3,N+1)
    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]
        z1 = result3[end]
        S+=s
        x_new = x1 + s*f1(x1,y1,z1,S)
        push!(result1,x_new)

        y_new = y1 + s*f2(x1,y1,z1,S)
        push!(result2,y_new)

        z_new = z1 + s*f3(x1,y1,z1,S)
        push!(result3,z_new)
    end
    return result1, result2, result3
end

#euler_4も作らなければ
function euler_4(f1,f2,f3,f4,x0,y0,z0,w0,stepsize,T,N)
    s = T(stepsize)
    result1 = T[x0]
    result2 = T[y0]
    result3 = T[z0] 
    result4 = T[w0]
    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    sizehint!(result3,N+1)
    sizehint!(result4,N+1)
    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]
        z1 = result3[end]
        w1 = result4[end]
        S+=s
        x_new = x1 + s*f1(x1,y1,z1,w1,S)
        push!(result1,x_new)

        y_new = y1 + s*f2(x1,y1,z1,w1,S)
        push!(result2,y_new)

        z_new = z1 + s*f3(x1,y1,z1,w1,S)
        push!(result3,z_new)

        w_new = w1 + s*f4(x1,y1,z1,w1,S)
        push!(result4,w_new)
    end
    return result1, result2, result3, result4
end

#ルンゲ・クッタ法
function rk4_1(f1,x0,stepsize,T,N)
    s = T(stepsize)
    result1 = T[x0]
    sizehint!(result1,N+1)
    S= T(0.0)

    for i in range(1,N)
        x1 = result1[end]
        S+=s
        k1 = f1(x1,S)
        k2 = f1(x1+k1*(s/2),S+(s/2))
        k3 = f1(x1+k2*(s/2),S+(s/2))
        k4 = f1(x1+(k3*s),S+s)
        x_new = x1 + (s)*(k1/6 + k2/3 + k3/3 + k4/6)
        push!(result1,x_new)
    end
    return result1
end

function rk4_2(f1,f2,x0,y0,stepsize,T,N)
    s = T(stepsize)

    result1 = [x0]
    result2 = [y0]

    sizehint!(result1,N+1)
    sizehint!(result2,N+1)

    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]

        S+=s

        k1 = T(f1(x1,y1,S))
        j1 = T(f2(x1,y1,S))

        k2 = T(f1(x1+k1*(s/2),y1+j1*(s/2),S+(s/2)))
        j2 = T(f2(x1+k1*(s/2),y1+j1*(s/2),S+(s/2)))

        k3 = T(f1(x1+k2*(s/2),y1+j2*(s/2),S+(s/2)))
        j3 = T(f2(x1+k2*(s/2),y1+j2*(s/2),S+(s/2)))

        k4 = T(f1(x1+k3*s,y1+j3*s,S+s))
        j4 = T(f2(x1+k3*s,y1+j3*s,S+s))

        x_new = x1 + (s)*(k1/6 + k2/3 + k3/3 + k4/6)
        push!(result1,x_new)

        y_new = y1 + (s)*(j1/6 + j2/3 + j3/3 + j4/6)
        push!(result2,y_new)
    end
    return result1, result2
end

function rk4_3(f1,f2,f3,x0,y0,z0,stepsize,T,N)
    s = T(stepsize)

    result1 = T[x0]
    result2 = T[y0]
    result3 = T[z0]

    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    sizehint!(result3,N+1)

    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]
        z1 = result3[end]

        S+=s

        k1 = f1(x1,y1,z1,S)
        j1 = f2(x1,y1,z1,S)
        l1 = f3(x1,y1,z1,S)

        k2 = f1(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),S+(s/2))
        j2 = f2(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),S+(s/2))
        l2 = f3(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),S+(s/2))

        k3 = f1(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),S+(s/2))
        j3 = f2(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),S+(s/2))
        l3 = f3(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),S+(s/2))

        k4 = f1(x1+k3*s,y1+j3*s,z1+l3*s,S+s)
        j4 = f2(x1+k3*s,y1+j3*s,z1+l3*s,S+s)
        l4 = f3(x1+k3*s,y1+k3*s,z1+l3*s,S+s)

        x_new = x1 + (s)*(k1/6 + k2/3+ k3/3 + k4/6)
        push!(result1,x_new)

        y_new = y1 + (s)*(j1/6 +j2/3 + j3/3 + j4/6)
        push!(result2,y_new)

        z_new = z1 + (s)*(l1/6 + l2/3 + l3/3 + l4/6)
        push!(result3,z_new)
    end
    return result1, result2, result3
end

#rk4_4も作らなければ

function rk4_4(f1,f2,f3,f4,x0,y0,z0,w0,stepsize,T,N)
    s = T(stepsize)
    result1 = T[x0]       
    result2 = T[y0]
    result3 = T[z0]
    result4 = T[w0]

    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    sizehint!(result3,N+1)
    sizehint!(result4,N+1)

    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]
        z1 = result3[end]
        w1 = result4[end]
        S+=s

        k1 = f1(x1,y1,z1,w1,S)
        j1 = f2(x1,y1,z1,w1,S)
        l1 = f3(x1,y1,z1,w1,S)
        m1 = f4(x1,y1,z1,w1,S)

        k2 = f1(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2))
        j2 = f2(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2))
        l2 = f3(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2))
        m2 = f4(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2))

        k3 = f1(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2))
        j3 = f2(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2))
        l3 = f3(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2))
        m3 = f4(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2))

        k4 = f1(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s)
        j4 = f2(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s)
        l4 = f3(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s)
        m4 = f4(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s)

        x_new = x1 + (s)*(k1/6 + k2/3+ k3/3 + k4/6)
        push!(result1,x_new)

        y_new = y1 + (s)*(j1/6 + j2/3 + j3/3 + j4/6)
        push!(result2,y_new)

        z_new = z1 + (s)*(l1/6 + l2/3 + l3/3 + l4/6)
        push!(result3,z_new)

        w_new = w1 + (s)*(m1/6 + m2/3 + m3/3 + m4/6)
        push!(result4,w_new)
    end
    return result1, result2, result3, result4
end


