
using LCIO
using Plots
using Glob
using LsqFit
using StatsBase
using Distributions
gr();

# read the files for calibration
fileList = readdir(glob"K0L*barrel*", "SingleParticles")
# build the list of input energies
energyString = r"_(\d+)GeV"
# get (match) the energy string from the filename, convert (parse) from string to Int16, build a dictionary from ints to filenames
energyMap = Dict((parse(Int16, match(energyString, fn)[1]), fn) for fn in fileList)
for (energy, filename) in energyMap
    println(energy, " GeV:\t", filename)
end

function getHitsFromFile(filename)
    eCalEnergies = Float64[]
    hCalEnergies = Float64[]
    eCalLengths = Int16[]
    hCalLengths = Int16[]
    LCIO.open(filename) do reader
        for (idx, event) in enumerate(reader)
            if idx > 10000
                break
            end
            hCalEnergy = 0.0
            hCalLength = 0
            eCalEnergy = 0.0
            eCalLength = 0
            # sum up the uncalibrated HCalHits
            for hit in getCollection(event, "HcalBarrelHits")
                hCalEnergy += getEnergy(hit)
                hCalLength += 1
            end
            # sum up the uncalibrated ECalHits
            # this needs to be sorted by layer, so we need a decoder
            EcalBarrelHits = getCollection(event, "EcalBarrelHits")
            decode = CellIDDecoder(EcalBarrelHits)
            for hit in EcalBarrelHits
                # calibrate the hits in the later layers with a higher number,
                # because they are behind thicker tungsten slabs
                factor = decode(hit)["layer"] < 20 ? 1 : 2
                eCalEnergy += factor*getEnergy(hit)
                eCalLength += 1
            end
            # fixme: Simple outlier cut
            if eCalLength+hCalLength < 100
                continue
            end
            push!(eCalEnergies, eCalEnergy)
            push!(eCalLengths, eCalLength)
            push!(hCalEnergies, hCalEnergy)
            push!(hCalLengths, hCalLength)
        end
    end
    return eCalEnergies, eCalLengths, hCalEnergies, hCalLengths
end;

eHits  = Dict{Int16, Vector{Float64}}()
eCount = Dict{Int16, Vector{Int16}}()
hHits  = Dict{Int16, Vector{Float64}}()
hCount = Dict{Int16, Vector{Int16}}()
for (energy, filename) in energyMap
    if energy < 10
        continue
    end
    println("Processing file for ", energy, " GeV")
    eCal, nEhits, hCal, nHhits = getHitsFromFile(filename)
    eHits[energy] = eCal
    eCount[energy] = nEhits
    hHits[energy] = hCal
    hCount[energy] = nHhits
end

# removes the 10% of the furthest outliers on either side
# no assumption about smoothness
function removeTails(distribution, cutOff=10)
    sort!(distribution)
    l = length(distribution)
    lcut = round(Int64, l * cutOff/100)
    hcut = round(Int64, l * (100-cutOff)/100)
    # start out with the whole distribution
    minDist = distribution[end] - distribution[1]
    low = 1
    high = l
    for idx = 1:lcut
        dist = distribution[hcut+idx] - distribution[idx]
        if dist < minDist
            minDist = dist
            low = idx+1
            high = hcut+idx
        end
    end
    return low, high
end;

function fitter(ecal, hcal, energy)
    function model(x, p)
        return p[1] .* x[1] + p[2] .* x[2]
    end
    fit = curve_fit(model, [ecal, hcal], energy, [50.0, 35.0])
    return fit.param
end
hCalCalibration = Dict{Int16, Float64}()
eCalCalibration = Dict{Int16, Float64}()
for energy in keys(hHits)
    ecal = eHits[energy]
    hcal = hHits[energy]
    calibration = fitter(ecal, hcal, energy)
    eCalCalibration[energy] = calibration[1]
    hCalCalibration[energy] = calibration[2]
    println(energy, " GeV: \t", eCalCalibration[energy], "\t", hCalCalibration[energy])
end

histogram([eCalCalibration[energy] .* eHits[energy] .+ hCalCalibration[energy] .* hHits[energy] for energy in keys(hHits)], fillalpha=0.5, linewidth=0, label=map(string, keys(hHits)))

# this function attempts a global fit and minimizes the offset b of the linear form y=mx+b
# parameters are ecal energies (×2 for the hits in the outer layers), hcal energies, particle energy
function lineFitter(ecal, hcal, truValues)
    function model(x, p)
        energies = Dict{Int16, Float64}()
        for e in truValues
            calibrated = p[1] .* x[1][e] + p[2] .* x[2][e]
            # cut the tails, fit a Normal distribution to the result
            low, high = removeTails(calibrated)
            n = Distributions.fit(Normal, calibrated[low:high])
            energies[e] = n.μ
        end
        # we are optimizing for the ratio of reconstructed energies to true values
        return [energies[e]/e for e in truValues]
    end
    fit = curve_fit(model, [ecal, hcal], 1.0, [50.0, 36.0])
    return fit.param
end
ECal, HCal = lineFitter(eHits, hHits, keys(hHits))
println("ECal calibration constant: ", ECal)
println("ECal calibration constant: ", HCal)

histogram([ECal .* eHits[energy] .+ HCal .* hHits[energy] for energy in keys(hHits)], fillalpha=0.5, linewidth=0, label=map(string, keys(hHits)))

y = Vector{Float64}()
yerr = Vector{Float64}()
x = Vector{Float64}()
for e in keys(hHits)
    d = ECal .* eHits[e] .+ HCal .* hHits[e]
    low, high = removeTails(d)
    gauss = Distributions.fit(Normal, d[low:high])
    push!(y, gauss.μ)
    push!(yerr, gauss.σ)
    push!(x, 1.0*e)
end
line(x, p) = p[1] + p[2]*x
l = curve_fit(line, x, y, [0.0, 1.0])
plot(x, y, yerr=yerr, marker=stroke(2), line=false, label="data")
leg = @sprintf("y=%.2f + %.2f*x", l.param[1], l.param[2])
plot!(x, l.param[1]+l.param[2]*x, label=leg)


