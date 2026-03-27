#using Pkg;Pkg.add("Plots")
using Plots
using Quadmath
n=5
# Define the input text
input_text = """
dot(p)_11(L=4) = 1.63900568901606809186556613225544795e-02
dot(p)_11(L=5) = 1.30268178079060740675110196017887465e-02
dot(p)_11(L=6) = 1.17260655278791101967879822430823147e-02
dot(p)_11(L=7) = 1.11901601319418966621021445536923404e-02
dot(p)_11(L=8) = 1.09529341304288466093308155978025843e-02
dot(p)_11(L=9) = 1.08401543967000342855719948633546335e-02
dot(p)_11(L=10) = 1.0784540602185199088126676050750971e-02
dot(p)_11(L=11) = 1.07585528937559276723828744641195542e-02
dot(p)_11(L=12) = 1.07497317313297676886314083755912569e-02
dot(p)_11(L=13) = 1.07516586171052854217543229090819615e-02
dot(p)_11(L=14) = 1.07605923569272988373972373854364713e-02
dot(p)_11(L=15) = 1.07741496176514345546366540760971751e-02
dot(p)_11(L=16) = 1.07907312095680612429580759595390448e-02
dot(p)_11(L=17) = 1.08092352143498924034516116161167956e-02
dot(p)_11(L=18) = 1.08288928270343397341581062418057192e-02
dot(p)_11(L=19) = 1.08491649146967476481099010862935376e-02
dot(p)_11(L=20) = 1.08696731238654590665924499540228919e-02
dot(p)_11(L=21) = 1.08901526719030887372437474596779213e-02
dot(p)_11(L=22) = 1.0910419511248370786893282691472279e-02
dot(p)_11(L=23) = 1.0930347279845586814175046298932498e-02
dot(p)_11(L=24) = 1.09498510048009874893167667363124184e-02
dot(p)_11(L=25) = 1.09688755043936413715531749681534663e-02
dot(p)_11(L=26) = 1.09873870805666424816396425168562127e-02
dot(p)_11(L=27) = 1.10053675310330081952979435773110025e-02
dot(p)_11(L=28) = 1.10228098077934538437781275281246568e-02
dot(p)_11(L=29) = 1.10397148525124689747028429013784921e-02
dot(p)_11(L=30) = 1.10560892790943219622375093434520577e-02
dot(p)_11(L=31) = 1.10719436703722571817196371283469021e-02
dot(p)_11(L=32) = 1.1087291322884949926475348728968633e-02
dot(p)_11(L=33) = 1.11021473205978289374030735926396412e-02
dot(p)_11(L=34) = 1.11165278514437880946348140527866322e-02
dot(p)_11(L=35) = 1.11304497039853683633651838171978831e-02
dot(p)_11(L=36) = 1.11439298982501065882840319856858821e-02
dot(p)_11(L=37) = 1.11569854168552853346190803163646117e-02
dot(p)_11(L=38) = 1.11696330112908990128785873089057353e-02
dot(p)_11(L=39) = 1.1181889064623183854830083759317368e-02
dot(p)_11(L=40) = 1.11937694965817866185748953238848942e-02
dot(p)_11(L=41) = 1.12052897004709130430517492597775508e-02
dot(p)_11(L=42) = 1.12164645039317176586085144668056291e-02
dot(p)_11(L=43) = 1.12273081475179325682301817740608205e-02
dot(p)_11(L=44) = 1.12378342765007651341946627651167682e-02
dot(p)_11(L=45) = 1.12480559424167345904177906214804747e-02
dot(p)_11(L=46) = 1.12579856117040698501604436945964788e-02
dot(p)_11(L=47) = 1.1267635179406149127684913054023582e-02
dot(p)_11(L=48) = 1.12770159864034073532616596321230524e-02
dot(p)_11(L=49) = 1.12861388390046940079679671938676995e-02
dot(p)_11(L=50) = 1.12950140300125003300421023113446041e-02
dot(p)_11(L=51) = 1.13036513605942507245153845424931907e-02
dot(p)_11(L=52) = 1.13120601624593923092160674093289078e-02
dot(p)_11(L=53) = 1.13202493199709974835775842408087002e-02
dot(p)_11(L=54) = 1.13282272919198968850139674783776542e-02
dot(p)_11(L=55) = 1.13360021327657618594894673992233752e-02
dot(p)_11(L=56) = 1.13435815132082487324574337541118043e-02
dot(p)_11(L=57) = 1.13509727399962924998503860026728501e-02
dot(p)_11(L=58) = 1.13581827749179719554542514645716878e-02
dot(p)_11(L=59) = 1.13652182529394407255063834282472021e-02
dot(p)_11(L=60) = 1.13720854994810733030021875063140864e-02
dot(p)_11(L=61) = 1.13787905468336447225469298060006896e-02
dot(p)_11(L=62) = 1.13853391497281623943186952438706447e-02
dot(p)_11(L=63) = 1.13917368000807686536164549934861366e-02
dot(p)_11(L=64) = 1.13979887409396120604238714390246946e-02
"""

#Extrapolation
const delta = 1
const Lmin = 4

# Extract the numbers using a regular expression
f = [parse(Float128, match(r"[-+]?\d*\.\d+([eE][-+]?\d+)?", line).match) for line in split(input_text, '\n') if match(r"[-+]?\d*\.\d+([eE][-+]?\d+)?", line) !== nothing]

#print(f)


function R0(f)
    R0f = Array{Float128}(undef, length(f) - delta)
    for i in 1:length(f)-delta
        L = Lmin - 1 + i
        R0f[i] = Float128(L / delta) * (f[i+delta] - f[i])
    end
    return R0f
end

function Rnu(nu, f)
    Rnuf = Array{Float128}(undef, length(f) - delta)
    for i in 1:length(f)-delta
        L = Lmin - 1 + i
        Rnuf[i] = f[i] + Float128(L / (nu * delta)) * (f[i+delta] - f[i])
    end
    return Rnuf
end



function extrapolationf(f)
    final = Rnu(1, Rnu(1, f))
    for l in 1:length(final)
        @inbounds begin
            L = l - 1 + Lmin
            final[l] = L * (f[l] - final[l])

        end
    end
    return final
end

# scale by 1/L
for i in 1:length(f)
    L = Lmin - 1 + i
    f[i] = f[i] / L
end
println(" -------------------------------------------------------- ")
final1 = extrapolationf(f)
function perform_extrapolation(f,n)
    result = f
    for i in 1:n
        result = Rnu(i, result)
        result = Rnu(i, result)
    end
    return result
end


final = perform_extrapolation(final1,n)
L_values = []
f_values = []

for l in 1:length(final)
    L = l - 1 + Lmin    
    if L > 40
        println(' ')
        println("Final(L=$L)= ", final[l])
        push!(L_values, L)
        push!(f_values, final[l])
    end
end
println(" -------------------------------------------------------- ")
println("Extrapolation up to R_$n= ",f_values[end])
# 1/L
inv_L = 1.0 ./ L_values

# 
plot(inv_L, f_values,
    seriestype=:scatter,
    xlabel="a/L",
    ylabel="d r1 / dz",
    label="Data Points",
    markersize=2,
    markercolor=:black,
    title="d r1 / dz vs. a/L",
    legend=:bottomleft,
    grid=true,
    dpi=300)
xlims!(minimum(inv_L) * 0.95, 1.05 * maximum(inv_L))
ylims!(0.9999*minimum(f_values), 1.00001*maximum(f_values))
baseval = 0.01200
hline!([baseval], label="baseline=$baseval", linestyle=:dash, color=:red)

#slope = (log(f_values[end]) - log(f_values[end-1])) / (inv_L[end] - inv_L[end-1])
#yint = slope*(0-inv_L[end]) + log(f_values[end])
#println("Extrapolation= ",exp(yint))
# 保存为高分辨率图片（可选）
savefig("dr1v1.png")

# 显示图形
display(plot!())
