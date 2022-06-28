function process_data()
	vars = matread("$path/data/DataToRun_pop09_boot.mat")
	data12 = vars["data12nopop"]
	data09 = vars["data09nopop"]
	bp     = vars["bpnopop"]

	# Fix type of dictionary for code performance
	data = Dict{Symbol,Dict{Symbol,Union{Vector{<:Number},Int64}}}()

	for (d_k, dataset) in zip([:of_09, :on_09, :on_12], [bp, data09, data12])
		data[d_k] = Dict{Symbol,Union{Vector{<:Number},Int64}}()
		# Type cast integers
		for k in ["cdindex","first","cdid","mktsize","popular"]
	           data[d_k][Symbol(k)] = vecI64(dataset[k])
	    end
		# Type cast floats
		for k in ["numlist","localint","conditiondif","p","obsweight","condition"]
	           data[d_k][Symbol(k)] = vecF64(dataset[k])
	   	end
		data[d_k][:N] = Int(dataset["N"])
		data[d_k][:M] = Int(dataset["M"])

		if d_k != :of_09
			data[d_k][:basecond] = vecF64(dataset["basecond"])
			data[d_k][:pdif]     = vecF64(dataset["pdif"])
		end
		if d_k == :on_12
			data[d_k][:disappear] = vecF64(dataset["disappear"])
		end
	end
	@save "$path/../data/input/data_to_run.jld2" data
	return data
end