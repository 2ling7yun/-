using Statistics
using Optim
#定义多项式函数
function polynomial(x::AbstractVector,ω::AbstractVector)
    M = length(ω) - 1;
    n = length(x); #x多长就有几个观测值
    X = ones(n) #n个元素的全1向量
    for i in 1:M
        X = hcat(X,x.^i);#按行并排
    end
    ŷ = X * ω;
    return ŷ
end
#定义误差方程
function error(y::AbstractVector,x::AbstractVector,ω::AbstractVector)
    ϵ = (y .-polynomial(x,ω)).^2; #误差向量求了平方
    ϵ = 1/2 * sum(ϵ) #求和然后除以2，即最小二乘法
    return ϵ
end

function polynomial_curve_fitting(y::AbstractVector,x::AbstractVector,M::Integer)
    ω0 = rand(M+1);#设置一个起始值，rand为生成一列随机数
    error_function(ω) = error(y,x,ω) #一个关于ω的函数
    result = optimize(error_function,ω0);
    ω̂ = Optim.minimizer(result);#找到最小值的过程
    return ω̂
end

# 测试
x = 0:0.05:1;
y = sin.(2*pi.*x) .+ rand(length(x))/20;
M = 3;
ω̂ = polynomial_curve_fitting(y,x,M)

#画图
using Plots
plot(x,y,seriestype=:scatter)
fullx = 0:0.01:1
plot!(fullx,polynomial(fullx,ω̂))
