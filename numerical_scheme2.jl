using FixedPointNumbers
using Statistics
#= 確率的丸めを使って数値計算を実装する =#

#確率的丸めを定義する関数を作成
function sr(x,digits,random)
    T = typeof(x)
    d = digits
    r = random
    e = eps(T)
    if (x - floor(x,digits=d)) / eps(T) <= r
        return T(floor(x,digits=d))
    else
        return T(floor(x,digits=d)) + e
    end
end

#オイラー法
function euler_1sr(f1,x0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type
    s = stepsize
    d = digsits
    result1 = [x0]
    S= T(0.0)
    r = rand(Uniform(0,1),itr) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)
    for i in range(1,itr)
        x1 = result1[end]
        S+=s
        x_new = x1 + sr(s*f1(x1,S),d,r[i])
        push!(result1,x_new)
    end
    return result1
end

function euler_2sr(f1,f2,x0,y0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type
    s = stepsize
    d = digsits
    result1 = [x0]
    result2 = [y0]
    S= T(0.0)
    r = rand(Uniform(0,1),itr,2) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)
    for i in range(1,itr)
        x1 = result1[end]
        y1 = result2[end]
        S+=s
        x_new = x1 + sr(s*f1(x1,y1,S),d,r[i,1])
        y_new = y1 + sr(s*f2(x1,y1,S),d,r[i,2])
        push!(result1,x_new)
        push!(result2,y_new)
    end
    return result1,result2
end

function euler_3sr(f1,f2,f3,x0,y0,z0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type    
    s = stepsize
    d = digsits
    result1 = [x0]
    result2 = [y0] 
    result3 = [z0]
    S= T(0.0)
    r = rand(Uniform(0,1),itr,3) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)
    for i in range(1,itr)
        x1 = result1[end]
        y1 = result2[end]
        z1 = result3[end]
        S+=s
        x_new = x1 + sr(s*f1(x1,y1,z1,S),d,r[i,1])
        y_new = y1 + sr(s*f2(x1,y1,z1,S),d,r[i,2])
        z_new = z1 + sr(s*f3(x1,y1,z1,S),d,r[i,3])
        push!(result1,x_new)
        push!(result2,y_new)
        push!(result3,z_new)
    end
    return result1,result2,result3
end

#euler_4sr

#ルンゲ・クッタ法
function rk4_1sr(f1,x0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type
    s = stepsize    
    d = digsits
    result1 = [x0]  
    S= T(0.0)   
    r = rand(Uniform(0,1),itr,4) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)

    for i in range(1,itr)
        x1 = result1[end]
        S+=s
        k1 = sr(s*f1(x1,s),d,r[i,1])
        k2 = sr(s*f1(x1+k1*(s/2),S+(s/2)),d,r[i,2])
        k3 = sr(s*f1(x1+k2*(s/2),S+(s/2)),d,r[i,3])
        k4 = sr(s*f1(x1+(k3*s),S+s),d,r[i,4])
        x_new = x1 + (k1/6 + k2/3+ k3/3 + k4/6)
        push!(result1,x_new)
    end

    return result1
end

function rk4_2sr(f1,f2,x0,y0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type
    s = stepsize
    d = digsits
    result1 = [x0]
    result2 = [y0]
    S= T(0.0)
    r = rand(Uniform(0,1),itr,8) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)
    
    for i in range(1,itr)
        x1 = result1[end]
        y1 = result2[end]
        S+=s
        k1 = sr(s*f1(x1,y1,S),d,r[i,1])
        j1 = sr(s*f2(x1,y1,S),d,r[i,2])

        k2 = sr(s*f1(x1+k1*(s/2),y1+j1*(s/2),S+(s/2)),d,r[i,3])
        j2 = sr(s*f2(x1+k1*(s/2),y1+j1*(s/2),S+(s/2)),d,r[i,4])

        k3 = sr(s*f1(x1+k2*(s/2),y1+j2*(s/2),S+(s/2)),d,r[i,5])
        j3 = sr(s*f2(x1+k2*(s/2),y1+j2*(s/2),S+(s/2)),d,r[i,6])

        k4 = sr(s*f1(x1+k3*s,y1+j3*s,S+s),d,r[i,7])
        j4 = sr(s*f2(x1+k3*s,y1+j3*s,S+s),d,r[i,8])

        x_new = x1 + (k1/6 + k2/3+ k3/3 + k4/6)
        push!(result1,x_new) 

        y_new = y1 + (j1/6 + j2/3+ j3/3 + j4/6) 
        push!(result2,y_new)
    end

    return result1,result2
end

function rk4_3sr(f1,f2,f3,x0,y0,z0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type
    s = stepsize
    d = digsits
    result1 = [x0]
    result2 = [y0]
    result3 = [z0]
    S= T(0.0)
    r = rand(Uniform(0,1),itr,12) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)
    for i in range(1,itr)
        x1 = result1[end]   
        y1 = result2[end]
        z1 = result3[end]
        S+=s
        k1 = sr(s*f1(x1,y1,z1,S),d,r[i,1])
        j1 = sr(s*f2(x1,y1,z1,S),d,r[i,2])  
        l1 = sr(s*f3(x1,y1,z1,S),d,r[i,3])

        k2 = sr(s*f1(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),S+(s/2)),d,r[i,4])   
        j2 = sr(s*f2(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),S+(s/2)),d,r[i,5])
        l2 = sr(s*f3(x1+k1/2,y1+(s/2),z1+(s/2),S+(s/2)),d,r[i,6])

        k3 = sr(s*f1(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),S+(s/2)),d,r[i,7])
        j3 = sr(s*f2(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),S+(s/2)),d,r[i,8])
        l3 = sr(s*f3(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),S+(s/2)),d,r[i,9])

        k4 = sr(s*f1(x1+k3*s,y1+j3*s,z1+l3*s,S+s),d,r[i,10])
        j4 = sr(s*f2(x1+k3*s,y1+j3*s,z1+l3*s,S+s),d,r[i,11])
        l4 = sr(s*f3(x1+k3*s,y1+j3*s,z1+l3*s,S+s),d,r[i,12])    

        x_new = x1 + (k1/6 + k2/3+ k3/3 + k4/6)
        push!(result1,x_new) 

        y_new = y1 + (j1/6 + j2/3+ j3/3 + j4/6) 
        push!(result2,y_new)

        z_new = z1 + (l1/6 + l2/3+ l3/3 + l4/6)
        push!(result3,z_new)
    end 

    return result1,result2,result3
end

#rk4_4sr
function rk_4sr(f1,f2,f3,f4,x0,y0,z0,w0,stepsize,Type,itr,digsits) #(関数，初期値，ステップサイズ，型(整数型など),反復回数，有効数字(10進数での))
    Random.seed!(1234) # 乱数のシードを指定する
    T = Type
    s = stepsize
    d = digsits
    result1 = [x0]
    result2 = [y0]
    result3 = [z0]
    result4 = [w0]
    S= T(0.0)
    r = rand(Uniform(0,1),itr,16) # 乱数の発生器を指定し乱数を発生させる(一様分布を用いているがベータ分布でも良い？)
    for i in range(1,itr)
        x1 = result1[end]
        y1 = result2[end]
        z1 = result3[end]
        w1 = result4[end]
        S+=s

        k1 = sr(s*f1(x1,y1,z1,w1,S),d,r[i,1])
        j1 = sr(s*f2(x1,y1,z1,w1,S),d,r[i,2])  
        l1 = sr(s*f3(x1,y1,z1,w1,S),d,r[i,3])
        m1 = sr(s*f4(x1,y1,z1,w1,S),d,r[i,4])

        k2 = sr(s*f1(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2)),d,r[i,5])   
        j2 = sr(s*f2(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2)),d,r[i,6])
        l2 = sr(s*f3(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2)),d,r[i,7])
        m2 = sr(s*f4(x1+k1*(s/2),y1+j1*(s/2),z1+l1*(s/2),w1+m1*(s/2),S+(s/2)),d,r[i,8])

        k3 = sr(s*f1(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2)),d,r[i,9])
        j3 = sr(s*f2(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2)),d,r[i,10])
        l3 = sr(s*f3(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2)),d,r[i,11])
        m3 = sr(s*f4(x1+k2*(s/2),y1+j2*(s/2),z1+l2*(s/2),w1+m2*(s/2),S+(s/2)),d,r[i,12])

        k4 = sr(s*f1(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s),d,r[i,13])
        j4 = sr(s*f2(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s),d,r[i,14])
        l4 = sr(s*f3(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s),d,r[i,15])
        m4 = sr(s*f4(x1+k3*s,y1+j3*s,z1+l3*s,w1+m3*s,S+s),d,r[i,16])    

        x_new = x1 + (k1/6 + k2/3+ k3/3 + k4/6)
        push!(result1,x_new) 

        y_new = y1 + (j1/6 + j2/3+ j3/3 + j4/6) 
        push!(result2,y_new)

        z_new = z1 + (l1/6 + l2/3+ l3/3 + l4/6)
        push!(result3,z_new)

        w_new = w1 + (m1/6 + m2/3+ m3/3 + m4/6)
        push!(result4,w_new)
    end
    
    return result1,result2,result3,result4
end