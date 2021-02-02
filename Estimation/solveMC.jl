using NLsolve

# set moving cost and time horizon vectors (first 12 elements are for type 1; last 12 elements are for type 2)
MC = -[-2.654795205, -3.144665335, -6.574142317, -7.064012447, -8.476504970, -8.966375099, -7.628074092, -8.117944222, -8.527064851, -9.016934981, -9.141404761, -9.631274891, -2.910429570, -3.400299700, -6.829776683, -7.319646812, -8.732139335, -9.222009465, -7.883708458, -8.373578588, -8.782699217, -9.272569347, -9.397039127, -9.886909257]
horzns = 65 .- [18, 18, 39, 39, 39, 39, 25, 25, 40, 40, 55, 55, 18, 18, 39, 39, 39, 39, 25, 25, 40, 40, 55, 55]
incs = [34530.18, 34530.18, 41364.76, 41364.76, 41364.76, 41364.76, 27142.27, 27142.27, 44071.81, 44071.81, 42667.49, 42667.49, 34530.18, 34530.18, 41364.76, 41364.76, 41364.76, 41364.76, 27142.27, 27142.27, 44071.81, 44071.81, 42667.49, 42667.49]

function executor(MCvec,horzns)
    # Present value formula
    function f!(F,x,a)
        # a = [beta, tau, T, MC]
        F[1] = sum( ( a[1]^(t-a[2]) )*log(x[1]) for t in a[2]:a[3] )-a[4]
    end

    beta=0.9
    tau=1
    for i=1:length(MCvec)
        g!(F,x) = f!(F,x,[beta, tau, horzns[i], MCvec[i]]) # create closure
        res = nlsolve(g!, [1.5])
        output = round.(res.zero[1]; digits=4)
        println("MC = ",MC[i],"; T = ",horzns[i],"; tau = ",tau,";wage ratio = ",output,";wage pct = ",round(100*(output-1);digits=1))
    end
end

executor(MC,horzns)

function executorPV(MCvec,horzns,incs)
    # Present value formula
    function f!(F,x,a)
        # a = [beta, tau, T, MC]
        F[1] = sum( ( a[1]^(t-a[2]) )*log(x[1]) for t in a[2]:a[3] )-a[4]
    end

    beta=0.9
    tau=1
    for i=1:length(MCvec)
        g!(F,x) = f!(F,x,[beta, tau, horzns[i], MCvec[i]]) # create closure
        res = nlsolve(g!, [1.5])
        output = round.(res.zero[1]; digits=4)
        println("MC = ",MC[i],"; T = ",horzns[i],"; tau = ",tau,";wage ratio = ",output,";PV = ",round( sum((beta^(t-1)) * ((res.zero[1]-1)*incs[i]) for t=1:horzns[i]); digits=0) )
    end
end

executorPV(MC,horzns,incs)

println("std dev of local amenities: ",sum((0.9^(t-1)) * (.077*incs[1]) for t=1:20) )
println("range of local amenities: ",sum((0.9^(t-1)) * (.302*incs[1]) for t=1:20) )
println("birth state bonus: ",sum((0.9^(t-1)) * (.189*incs[1]) for t=1:20) )
