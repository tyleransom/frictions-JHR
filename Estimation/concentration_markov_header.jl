
# Simulate an economy of F firms in J locations with varying levels of moving costs for workers to move firms
# Match the hypothetical parameter values to observed flows in the SIPP

@views @inbounds function datagenAmen(ew,sigw,siga,covwa,J=4,flag=true)
    Random.seed!(32)
    # this function generates wages and amenities for firms
    s0 = (1/J)*ones(J) # initial worker distribution
    Σ = [sigw covwa; covwa siga] # covariance between wages and amenities (should be neg. bec. compens. differentials)
    d = MvNormal(zeros(2),Σ)
    r = rand(d,J);
    η = round.(r[1,:]; digits=1)
    α = r[2,:]
    w = ew .+ η
    #w    = ew .+ ones(N)*ones(1,J)
    if flag==false
        α = 0 .* α
    end
    return s0,w,α,J
end

@views @inbounds function getθmkt(J,numMkt,th1,th2)
    #define markets
    @assert mod(J,numMkt)==0 "number of firms must be a multiple of the number of markets"
    firmpermkt = convert(Int64,J/numMkt)
    mktlist = []
    for m=1:numMkt
        append!(mktlist,m*ones(firmpermkt))
    end
    mktlist = convert(Array{Float64,1},mktlist)

    mktblkdiag = zeros(J,J)
    for m=1:numMkt
        mktblkdiag[findall(mktlist.==m),findall(mktlist.==m)] .= 1
    end

    # θ1: market moving cost
    θ1 = th1*(ones(J,J) - mktblkdiag)

    # θ2: firm switching cost
    θ2 = th2*(ones(J,J) - I)

    return θ1,θ2,mktlist,mktblkdiag
end

@views @inbounds function getdynshares(α,γ,θ1,θ2,w,J,β,T,ρ)
    #===
    α: J-vector of amenities
    γ: scalar for wage preference
    θ1: JxJ matrix of market moving costs
    θ2: JxJ matrix of firm switching costs
    w: J-vector of wages
    J: number of firms
    β: workers' discount factor
    T: time horizon
    ===#
    TransMat = zeros(J,J)
    if β==0
        for j=1:J # previous firm
            u = α .+ γ*w .+ θ1[j,:] .+ θ2[j,:]
            TransMat[j,:] .= exp.(u) ./ sum(exp.(u)) # Pr(k in t | j in t-1)
        end
    else
        v  = zeros(J)
        FV = zeros(T+1,J)
        for t=T:-1:1
            for k=1:J # previous firm
                v = α .+ γ*ρ^(t-1)*w .+ θ1[k,:] .+ θ2[k,:] .+ FV[t+1,:]
                if t==1
                    for j=1:J
                        v = α .+ γ*ρ^(t-1)*w .+ θ1[j,:] .+ θ2[j,:] .+ FV[t+1,:]
                        TransMat[j,:] .= exp.(v) ./ sum(exp.(v))
                    end
                end
                FV[t,k] = β*log(sum(exp.(v)))
            end
        end
    end
    return TransMat
end

@views @inbounds function get_εqw(step_size,TM0,α,γ,θ1,θ2,w,J,β,T,ρ)
    TM   = zeros(J,J)
    rec0 = zeros(J)
    rec1 = zeros(J)
    sep0 = zeros(J)
    sep1 = zeros(J)
    dw   = zeros(J)
    if β==0
        for k=1:J
            rec0[k] = sum(TM0[:,k])-TM0[k,k]
            sep0[k] = sum(TM0[k,:])-TM0[k,k]
            w1 = deepcopy(w)
            w1[k] += step_size
            dw[k] = w1[k]-w[k]
            for j=1:J
                u = α .+ γ*w1 .+ θ1[j,:] .+ θ2[j,:]
                TM[j,:] .= exp.(u) ./ sum(exp.(u))
            end
            rec1[k] = sum(TM[:,k])-TM[k,k]   
            sep1[k] = sum(TM[k,:])-TM[k,k]   
        end
    else
        v  = zeros(J)
        FV = zeros(T+1,J)
        for t=T:-1:1
            for k=1:J # previous firm
                v = α .+ γ*ρ^(t-1)*w .+ θ1[k,:] .+ θ2[k,:] .+ FV[t+1,:]
                if t==1
                    for m=1:J
                        rec0[m] = sum(TM0[:,m])-TM0[m,m]
                        sep0[m] = sum(TM0[m,:])-TM0[m,m]
                        w1 = deepcopy(w)
                        w1[m] += step_size
                        dw[m] = w1[m]-w[m]
                        for j=1:J
                            vv = α .+ γ*w1 .+ θ1[j,:] .+ θ2[j,:] .+ FV[t+1,:]
                            TM[j,:] .= exp.(vv) ./ sum(exp.(vv))
                        end
                        rec1[m] = sum(TM[:,m])-TM[m,m]   
                        sep1[m] = sum(TM[m,:])-TM[m,m]   
                    end
                end
                FV[t,k] = β*log(sum(exp.(v)))
            end
        end
    end
    εqw  = (sep1 .- sep0)./dw
    εrw  = (rec1 .- rec0)./dw
    return εqw,εrw
end

@views @inbounds function wrapper(exw,sgw,sga,cvwa,J,M,amenflag,wagegamma,movecost,switchcost,wagestep,discfac=0,timehorizon=40,ARcoef=.74)
    to = TimerOutput()
    @timeit to "elaps" begin
        s0,wage,amen,J = datagenAmen(exw,sgw,sga,cvwa,J,amenflag)
        θ_1,θ_2,MM,MM2 = getθmkt(J,M,movecost,switchcost) 
        Tmat = getdynshares(amen,wagegamma,θ_1,θ_2,wage,J,discfac,timehorizon,ARcoef)
        quit_elast,rec_elast = get_εqw(wagestep,Tmat,amen,wagegamma,θ_1,θ_2,wage,J,discfac,timehorizon,ARcoef) 
        wage_elast = rec_elast - quit_elast

        #migration rate
        migrat = sum(Tmat.*(1 .- MM2);dims=2) |> mean
        # should match .03411 in SIPP

        #job turnover rate
        switchrat = sum(Tmat.*(ones(size(Tmat)) - I);dims=2) |> mean
        # should match .2111054 in SIPP

        #get steady state of Markov Chain:
        capT = 20
        s = zeros(capT,J)
        for t=1:capT
            s[t,:]=s0'*Tmat^(t-1)
        end
    end
    result = DataFrame(NF = J, NM = M, avg_w = exw, sigma_w = sgw, sigma_amen = sga, cov_w_a = cvwa, gamma = convert(Float64,wagegamma), θ_move = convert(Float64,movecost), θ_switch = convert(Float64,switchcost), wagestep = wagestep, avgquiteps = mean(quit_elast), avgreceps = mean(rec_elast), avgelast = mean(wage_elast), medelast = median(wage_elast), stdelast = std(wage_elast), minelast = minimum(wage_elast), maxelast = maximum(wage_elast), migrate = migrat, switchrate = switchrat, avgshr = mean(s[capT,:]), minshr = minimum(s[capT,:]), maxshr = maximum(s[capT,:]), disc_fac = convert(Float64,discfac), comptime = TimerOutputs.tottime(to)/1e9, cor_w_elast = cor(wage,wage_elast), cor_a_elast = cor(amen,wage_elast))
    return result
end

function exporter(numfirm,nummkt,beta=0)
    if beta==0
        fname  = string("results_",nummkt,"_",numfirm,"_static.csv")
        fname2 = string("results_all_",nummkt,"_",numfirm,"_static.csv")
    else
        fname  = string("results_",nummkt,"_",numfirm,".csv")
        fname2 = string("results_all_",nummkt,"_",numfirm,".csv")
    end
    results = DataFrame()
    temp = wrapper(7.96,.1,.14,-.07,numfirm,nummkt,true,1,-50,-50,.0001,beta)
    append!(results,temp)
    CSV.write(fname, temp)
    for thmov in [-6:.5:0;]
        for thsw in [-6:.5:0;] 
            for gam = 1:5
                temp = wrapper(7.96,.1,.14,-.07,numfirm,nummkt,true,gam,thmov,thsw,.0001,beta)
                append!(results,temp)
                CSV.write(fname, temp; append=true)
            end
        end
    end
    CSV.write(fname2, results)
    return nothing
end

