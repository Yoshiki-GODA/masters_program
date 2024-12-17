using FixedPointNumbers,LinearAlgebra
include("numerical_scheme.jl")

function stormer_verlet(f1,f2,x0,y0,stepsize,T,N) #間違っている！！！
    s = T(stepsize)
    result1 = T[x0]
    result2 = T[y0]
    sizehint!(result1,N+1)
    S = T(0.0)
    for i in range(1,N)
        x1 = result1[end]
        y1 = result2[end]
        S+=s
        x_half = x1 + (s/2)*f1(x1,y1,S) 
        y_new = y1 + s*f2(x_half,y1,S)
        x_new = x_half + (s/2)*f1(x_half,y_new,S)
        push!(result1,x_new)
        push!(result2,y_new)
    end
    return result1, result2
end

function newton(H,dH,x_ini,stepsize,Type,itr,tol)
    T = Type
    x = x_ini
    
    for i in 1:itr
        if norm(H(x),2) < tol
            break
        end
        x = x - dH(x) \ H(x)
    end
    return x,norm(H(x))
end

function newton_shift(H,dH,x_ini,stepsize,shift,Type,itr,tol)
    T = Type
    x_shift = x_ini + shift
    
    for i in 1:itr
        if norm(H(x_shift),2) < tol
            break
        end
        x_shift = x_shift - dH(x_shift) \ H(x_shift)
    end
    x = x_shift
    return x,norm(H(x_shift),2)
end

function my_newton(H,dH_inv,x_ini,stepsize,Type,itr,tol)
    T = Type
    x = x_ini
    for i in 1:itr
        if norm(H(x),2) < tol
            break
        end
        x = x - dH_inv(x)*H(x)
    end
    return x, norm(H(x))
end

function spm_lv_T(v_0,f_1,f_2,vec_t,mat_inv_t,p,stepsize,Type,itr,itr_newton,tol_newton)
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
        v = newton(vec,mat,v_ini,s,T,itr_newton,tol_newton)[1]
        push!(result1,v[1])
        push!(result2,v[2])
    end
    return result1,result2
end

function spm_lv_T_shift(v_0,f_1,f_2,vec_t,mat_t,p,stepsize,shift,Type,itr,itr_newton,tol_newton)
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
        mat(v_new) = mat_t(v_new,[result1[end],result2[end]],s,p)
        v = newton_shift(vec,mat,v_ini,s,shift,T,itr_newton,tol_newton)[1]    
        push!(result1,v[1])
        push!(result2,v[2])
    end
    return result1,result2
end

function euler_rescale_1(f1,x0,stepsize,scaleratio,Type,itre)
    T = Type
    N = itre
    r = T(scaleratio)
    s = T(stepsize)
    result1 = T[x0]
    sizehint!(result1,N+1)
    S = T(0.0)
    for i in range(1,N)
        x1 = result1[end]/r
        S+=s
        x_new_scaled = x1 + (s/r)*f1(x1,S,r)
        x_new = r*x_new_scaled
        push!(result1,x_new)
    end
end

function euler_rescale_2(f1,f2,x0,y0,stepsize,scaleratio,Type,itre)
    T = Type
    N = itre
    r = T(scaleratio)
    s = T(stepsize)
    result1 = Type[x0]
    result2 = Type[y0]
    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]/r
        y1 = result2[end]/r
        S+=s
        x_new_scaled = x1 + (s*f1(x1,y1,S,r))/r
        x_new = r*x_new_scaled
        y_new_scaled = y1 + (s*f2(x1,y1,S,r))/r
        y_new = r*y_new_scaled
        push!(result1,x_new)
        push!(result2,y_new)
    end
    return result1, result2
end

function euler_rescale_3(f1,f2,f3,x0,y0,z0,stepsize,scaleratio,Type,itre)
    T = Type
    N = itre
    r = T(scaleratio)
    s = T(stepsize)
    result1 = T[x0]
    result2 = T[y0]
    result3 = T[z0]
    sizehint!(result1,N+1)
    sizehint!(result2,N+1)
    sizehint!(result3,N+1)
    S= T(0.0)
    for i in range(1,N)
        x1 = result1[end]/r
        y1 = result2[end]/r
        z1 = result3[end]/r
        S+=s
        x_new_scaled = x1 + (s*f1(x1,y1,z1,S,r))/r
        x_new = r*x_new_scaled
        push!(result1,x_new)

        y_new_scaled = y1 + (s*f2(x1,y1,z1,S,r))/r
        y_new = r*y_new_scaled
        push!(result2,y_new)

        z_new_scaled = z1 + (s*f3(x1,y1,z1,S,r))/r
        z_new = r*z_new_scaled
        push!(result3,z_new)
    end
    return result1, result2, result3
end
#=function vec_spm(v_new,v_old,s,p,Type)
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        x_new - x_old + s*((b*(exp(y_new)-exp(y_old))-a*(y_new - y_old))/(y_new - y_old)),
        y_new - y_old - s*((c*(exp(x_new)-exp(x_old))-d*(x_new - x_old))/(x_new - x_old))
    ]
end

function mat_spm(v_new,v_old,s,p)
    a,b,c,d = p
    x_old,y_old = v_old
    x_new,y_new = v_new
    return [
        1 s*b*((exp(y_new)/(y_new-y_old))-((exp(y_new)-exp(y_old))/((y_new-y_old)^2)))
        s*c*((exp(x_new)/(x_new-x_old))-((exp(x_new)-exp(x_old))/((x_new-x_old)^2))) 1
    ]
end


function H_lv(v_old,v_new,stepsize,Type)
    x_old,y_old = v_old
    x,y = v_new
    s = Type(stepsize)
    H = x + y -log(x) - log(y)
    return [
        (x_new - x_old)^2*(y_new - y_old)+s*(x_old)*(y_old)*(H(x_new,y_new)-H(x_new,y_old))
        (y_new - y_old)^2*(x_new - x_old)+s*(y_old)*(x_old)*(H(x_new,y_old)-H(x_old,y_old))
    ]
end

function dH_lv(v_old,v_new,stepsize,Type)
    x_old,y_old = v_old
    x,y = v_new
    s = Type(stepsize)
    return [
        2*(x-x_old)*(y-y_old)+s*(x_old)*(y_old)*((y-y_old)+log(y_old/y)) (x-x_old)^2+s*(x_old)*(y_old)*(T(1.0)-(T(1.0)/y))
        (y-y_old)^2+s*(x_old)*(y_old)*(y-y_old)*(T(1.0)-(T(1.0)/x)) 2*(y-y_old)*(x-x_old)+s*(y_old)*(x_old)*((x-x_old)+log(x_old/x))
    ]
end

function structure_preserving_lv(f1,f2,H,dH,x0,y0,Type,stepsize,itr,itr_n,tol)
    x,y = v #まだ途中
    T = Type
    s = T(stepsize)
    result1 = T[x0]
    result2 = T[y0]
    v_1 = newton(H, dH, [x0,y0], s, T, itr_n, tol)[1] 
    sizehint!(result1,itr+1)
    sizehint!(result2,itr+1)
    v1 = newton(H, dH, [x0,y0], s, T, itr_n, tol)[1]
    push!(result1,v1[1][1])
    push!(result2,v1[1][2])
    S = T(0.0)
    for i in range(itr-1)
        x1 = result1[end]
        y1 = result2[end]
        S+=s
        v = newton(H, dH, [x1,y1], s, T, itr_n, tol)[1]
        v_new = euler_2(f1,f2,v[1][1],v[1][2],s,T,1)
        push!(result1,v_new[1][1])
        push!(result2,v_new[2][2])
    end
    return result1,result2
end=#