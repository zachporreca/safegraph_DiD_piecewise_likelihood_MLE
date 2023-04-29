
#SIMULATING
using Pkg
Pkg.add("RData")
Pkg.add("CodecBzip2")
Pkg.add("Distributions")
Pkg.add("Optim")
Pkg.add("Distributed")
using RData
import CodecBzip2
rdata=load("") #filepath
rdata=rdata[""] #name of data.frame from R
v=rdata[:,1]
x1=rdata[:,2]
x2=rdata[:,3]
x3=rdata[:,4]

using Distributions, Optim



#Fixing for truncated optimization

function pwll(v, x1, x2, x3, a, b1, b2, b3)
    d = Poisson(exp(a+b1*x1+b2*x2+b3*x3))
    
    if v > 4 
        return pdf(d, v)
    elseif v == 4
        return pdf(d, 4) + pdf(d, 3) + pdf(d, 2)
    else
        return pdf(d, 1) + pdf(d, 0)
    end
end


fixILL(x) = iszero(x) ? 1e-17 : x

est(a, b1, b2, b3, v, x1, x2, x3) = -1 .* (sum(log10.(fixILL.(pwll.(v, x1, x2, x3, a, b1, b2, b3))))-cdf(Poisson(exp(a+b1*mean(x1)+b2*mean(x2)+b3*mean(x3))), 1))


obj(input; v = v, x1 = x1, x2 = x2, x3 = x3) = est(input[1], input[2], input[3], input[4], v, x1, x2, x3)

optimum=optimize(obj, [0, 0, 0, 0]) #starting values can be changed to increase speed
MLE=optimum.minimum
parameter=optimum.minimizer


#bootstrapping SE

intercepts_bs=Float64[]
post_bs=Float64[]
treat_bs=Float64[]
post_treat_bs=Float64[]

using Distributions, Optim, Distributed


@distributed for i in 1:100 #100 is number of bootstrap iterations here, 
    print(i)
    jsample=sample([1:1:29477729;],29477729, replace=true) #29477729 is number of observations
    jv=v[jsample]
    jx1=x1[jsample]
    jx2=x2[jsample]
    jx3=x2[jsample]



    #Fixing for truncated optimization
    
    function pwll(jv, jx1, jx2, jx3, a, b1, b2, b3)
        d = Poisson(exp(a+b1*jx1+b2*jx2+b3*jx3))
        
        if jv > 4 
            return pdf(d, jv)
        elseif jv == 4
            return pdf(d, 4) + pdf(d, 3) + pdf(d, 2)
        else
            return pdf(d, 1) + pdf(d, 0)
        end
    end
    
    
    fixILL(x) = iszero(x) ? 1e-17 : x
    
    est(a, b1, b2, b3, jv, jx1, jx2, jx3) = -1 .* (sum(log10.(fixILL.(pwll.(jv, jx1, jx2, jx3, a, b1, b2, b3))))-cdf(Poisson(exp(a+b1*mean(jx1)+b2*mean(jx2)+b3*mean(jx3))), 1))
    
    
    obj(input; jv = jv, jx1 = jx1, jx2 = jx2, jx3 = jx3) = est(input[1], input[2], input[3], input[4], jv, jx1, jx2, jx3)
    
    joptimum=optimize(obj, [0, 0, 0, 0]) #can change starting values
    jMLE=joptimum.minimum
    jparameter=joptimum.minimizer
    push!(intercepts_bs, jparameter[1])
    push!(post_bs, jparameter[2])
    push!(treat_bs, jparameter[3])
    push!(post_treat_bs, jparameter[4])

    print("percent")
end

using Statistics

std(intercepts_bs)
std(post_bs)
std(treat_bs)
std(post_treat_bs)

Pkg.add("DataFrames")
using DataFrames
bootstrapped = DataFrame(intercept=intercepts_bs, post=post_bs, treat=treat_bs, post_treat=post_treat_bs)

Pkg.add("CSV")
using CSV
CSV.write("", bootstrapped) #file path between quotation marks, saving a CSV of bootstrapped estimates
