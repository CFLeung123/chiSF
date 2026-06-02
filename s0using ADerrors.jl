using ADerrors
using Plots
using Quadmath      # provides Float128
using Printf

# ----------------------------------------------------------------------
# 1. Input data (high‑precision values parsed as Float128)
# ----------------------------------------------------------------------
input_text = """
p_11(L=4) = -4.49907573341427211096958142990506253e-02
p_11(L=5) = -4.65355795773412099869057322811882875e-02
p_11(L=6) = -4.81281731163087827556765268592643903e-02
p_11(L=7) = -4.95702422033655985577606656210849701e-02
p_11(L=8) = -5.08391340851986160495966132145604268e-02
p_11(L=9) = -5.19566845126150709412839952941867111e-02
p_11(L=10) = -5.29502877539142924220673610171327389e-02
p_11(L=11) = -5.38432019389609770522100207594702068e-02
p_11(L=12) = -5.46534766598364944139110703439455076e-02
p_11(L=13) = -5.53949054969890212736073722947356637e-02
p_11(L=14) = -5.6078133175487548732648101681164621e-02
p_11(L=15) = -5.67115173238096273257391488849784979e-02
p_11(L=16) = -5.73017379500439338154369826289987276e-02
p_11(L=17) = -5.78542208975699250978359405035678573e-02
p_11(L=18) = -5.83734352407383133901856404506284316e-02
p_11(L=19) = -5.88631064917523749218322788308746809e-02
p_11(L=20) = -5.93263728857670259806217939439489634e-02
p_11(L=21) = -5.97659023210444832634106252338485771e-02
p_11(L=22) = -6.01839814462311195162320275867457237e-02
p_11(L=23) = -6.05825845866426177874857440490802376e-02
p_11(L=24) = -6.09634277896220552443999258015373928e-02
p_11(L=25) = -6.13280117001846033954507348951347856e-02
p_11(L=26) = -6.16776559303643203993707780380306902e-02
p_11(L=27) = -6.20135268683083941294980368040016331e-02
p_11(L=28) = -6.23366603711289280353218905864791311e-02
p_11(L=29) = -6.26479804274725700160054456273773526e-02
p_11(L=30) = -6.29483146162200220445676724589323786e-02
p_11(L=31) = -6.32384069969094672990028748236828601e-02
p_11(L=32) = -6.35189289254388510190427727767853453e-02
p_11(L=33) = -6.37904881816949744700519762077691368e-02
p_11(L=34) = -6.40536367144822668475143148392480583e-02
p_11(L=35) = -6.43088772467657483707881183335058349e-02
p_11(L=36) = -6.45566689359908352852982773642949288e-02
p_11(L=37) = -6.47974322466104608259135865632370952e-02
p_11(L=38) = -6.50315531623828800693187816195007926e-02
p_11(L=39) = -6.52593868426115089708451269757736458e-02
p_11(L=40) = -6.54812608078703228448886512996917417e-02
p_11(L=41) = -6.56974777258326114749511993887405441e-02
p_11(L=42) = -6.59083178557914507409958188760142132e-02
p_11(L=43) = -6.61140412007112278510551963385682671e-02
p_11(L=44) = -6.63148894077069543296036369300421442e-02
p_11(L=45) = -6.65110874513447114271890125747750237e-02
p_11(L=46) = -6.67028451288061376838654479591252621e-02
p_11(L=47) = -6.68903583915379030899941979883610557e-02
p_11(L=48) = -6.7073810534336577959127608423223351e-02
p_11(L=49) = -6.72533732597598175720843257617582612e-02
p_11(L=50) = -6.7429207633194459705794990361505877e-02
p_11(L=51) = -6.76014649417613539447028086051833257e-02
p_11(L=52) = -6.77702874684233882877480941042133738e-02
p_11(L=53) = -6.79358091911289460423684309647149973e-02
p_11(L=54) = -6.80981564155205748530494470734811596e-02
p_11(L=55) = -6.82574483486293821466915081901484844e-02
p_11(L=56) = -6.8413797620027990375441450494623555e-02
p_11(L=57) = -6.85673107561028252898727631620868187e-02
p_11(L=58) = -6.87180886124086532032306287371694591e-02
p_11(L=59) = -6.88662267684668989982365527931737458e-02
p_11(L=60) = -6.90118158888496379285193340283563851e-02
p_11(L=61) = -6.91549420539409866080523450309942697e-02
p_11(L=62) = -6.92956870633766590624828854829992394e-02
p_11(L=63) = -6.94341287148220946745729804037483093e-02
p_11(L=64) = -6.957034106045256522839409528463877e-02
"""

lines = split(input_text, '\n')
f_big = [parse(Float128, match(r"[-+]?\d*\.\d+([eE][-+]?\d+)?", line).match)
         for line in lines if match(r"[-+]?\d*\.\d+([eE][-+]?\d+)?", line) !== nothing]

const delta = 1
const Lmin = 4
const Lmax = 64
const n_extrap = 5

# Error model: returns Float64 errors
function compute_errors(Lmin, Lmax, log10_error_min; slope=0.0539)
    Lmin_f64 = Float64(Lmin)
    slope_f64 = Float64(slope)
    log10_min = Float64(log10_error_min)
    b = log10_min - slope_f64 * Lmin_f64
    [10.0^(slope_f64 * L + b) for L in Lmin:Lmax]
end
errors_f64 = compute_errors(Lmin, Lmax, -29)

# ----------------------------------------------------------------------
# 2. High‑precision deterministic calculation (central values)
#    using Float128
# ----------------------------------------------------------------------
function R0_big(f::Vector{Float128})
    n = length(f) - delta
    result = Vector{Float128}(undef, n)
    for i in 1:n
        L = i - 1 + Lmin
        result[i] = (Float128(L) / delta) * (f[i+delta] - f[i])
    end
    return result
end

function Rnu_big(nu, f::Vector{Float128})
    n = length(f) - delta
    result = Vector{Float128}(undef, n)
    for i in 1:n
        L = i - 1 + Lmin
        result[i] = f[i] + (Float128(L) / (nu * delta)) * (f[i+delta] - f[i])
    end
    return result
end

function perform_extrapolation_big(f::Vector{Float128}, nsteps::Int)
    result = R0_big(f)
    for step in 1:nsteps
        result = Rnu_big(step, result)
        result = Rnu_big(step, result)
    end
    return result
end

final_central = perform_extrapolation_big(f_big, n_extrap)

# ----------------------------------------------------------------------
# 3. Error propagation using ADerrors.jl (Float64)
# ----------------------------------------------------------------------
f_central_f64 = Float64.(f_big)
uw_data = [uwreal([f_central_f64[i], errors_f64[i]], "L=$(Lmin+i-1)")
           for i in 1:length(f_big)]

function R0_err(f::Vector{uwreal})
    n = length(f) - delta
    res = Vector{uwreal}(undef, n)
    for i in 1:n
        L = i - 1 + Lmin
        res[i] = (L / delta) * (f[i+delta] - f[i])
    end
    return res
end

function Rnu_err(nu, f::Vector{uwreal})
    n = length(f) - delta
    res = Vector{uwreal}(undef, n)
    for i in 1:n
        L = i - 1 + Lmin
        res[i] = f[i] + (L / (nu * delta)) * (f[i+delta] - f[i])
    end
    return res
end

function perform_extrapolation_err(f::Vector{uwreal}, nsteps::Int)
    result = R0_err(f)
    for step in 1:nsteps
        result = Rnu_err(step, result)
        result = Rnu_err(step, result)
    end
    return result
end

final_err = perform_extrapolation_err(uw_data, n_extrap)

# ----------------------------------------------------------------------
# 4. Extract results for L > 40
# ----------------------------------------------------------------------
L_vals = Int[]
final_vals = Float128[]
final_errs = Float64[]

for i in 1:length(final_err)
    L = i - 1 + Lmin
    if L > 4+length(final_err)-6
        ADerrors.uwerr(final_err[i])
        push!(L_vals, L)
        push!(final_vals, final_central[i])
        push!(final_errs, ADerrors.err(final_err[i]))
    end
end

println("--------------------------------------------------------")
println("Extrapolation up to R_$n_extrap = ")
for i in eachindex(L_vals)
    @printf("L=%d: %.60f ± %.30f\n", L_vals[i], Float64(final_vals[i]), final_errs[i])
end
println("--------------------------------------------------------")

# ----------------------------------------------------------------------
# 5. Plot with clearly visible error bars
# ----------------------------------------------------------------------
inv_L = 1.0 ./ L_vals
final_vals_f64 = Float64.(final_vals)

# Use errorbar series type for better control
plot(inv_L, final_vals_f64,
     yerr = final_errs,        # 直接使用 yerr 关键字
     seriestype = :scatter,    # 系列类型保持为 scatter
     markersize = 2,
     markercolor = :black,
     xlabel = "a/L",
     ylabel = "s0",
     label = "Data Points",
     title = "s0 vs. a/L (high‑precision central values)",
     legend = :bottomleft,
     grid = true,
     dpi = 300)

hline!([-0.0084434319], label = "baseline=", linestyle = :dash, color = :red)
xlims!(minimum(inv_L)*0.95, maximum(inv_L)*1.05)
#ylims!(minimum(final_vals_f64)*1.00001, maximum(final_vals_f64)*0.999999)
savefig("s0SWv2_bigfloat.png")
display(plot!())