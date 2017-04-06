using LCIO
using Plots
using Glob
using LsqFit
gr()


function calibrateHCal(ecal, hcal, simulated)
	# for now, we'll only adjust the hcal calibration
	model(x, p) = ecal .+ p[1] .* x
	fit = curve_fit(model, hcal, simulated, [0.5])
	println(fit.param)
end

# read the files for calibration
fileList = readdir(glob"K0L*barrel*", "SingleParticles")
println(fileList)

# build the list of input energies
energyString = r"_(\d+)GeV"
# get the energy string from the filename, convert from string to int, build a dictionary from ints to filenames
energyMap = Dict((parse(Int64, match(energyString, fn)[1]), fn) for fn in fileList)
println(energyMap)

for (energy, fn) in energyMap
	singleParticleECal = Float64[]
	singleParticleHCal = Float64[]
	LCIO.open(fn) do reader
		for (idx, event) in enumerate(reader)
			if idx > 10000
				break
			end
			hCalEnergy = 0.0
			eCalEnergy = 0.0
			# for cn in getCollectionNames(event)
			# 	println(cn)
			# end
			# break
			# sum up the uncalibrated HCalHits
			for hit in getCollection(event, "HcalBarrelHits")
				hCalEnergy += getEnergy(hit)
			end
			# sum up the calibrated ECalHits
			for hit in getCollection(event, "EM_BARREL")
				eCalEnergy += getEnergy(hit)
			end
			push!(singleParticleHCal, hCalEnergy)
			push!(singleParticleECal, eCalEnergy)
		end
	end
	histogram(singleParticleHCal)
	gui()
	scatter(singleParticleECal, singleParticleHCal)
	gui()
	# for now only look at one energy
	calibrateHCal(singleParticleECal, singleParticleHCal, energy)
	break
end
