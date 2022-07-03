"""
```
function process_data()
```
Exists to convert + fix data types in the .mat file to a .jld2 file.
"""
function process_data()
	vars = matread("$INPUT/DataToRun_pop09_boot.mat")
	data12 = vars["data12nopop"]
	data09 = vars["data09nopop"]
	bp     = vars["bpnopop"]

	# Fix type of dictionary for code performance
	data = Dict{Symbol,Dict{Symbol,Vector{<:Number}}}()

	for (d_k, dataset) in zip([:of_09, :on_09, :on_12], [bp, data09, data12])
		data[d_k] = Dict{Symbol,Vector{<:Number}}()
		# Type cast integers
		for k in ["cdindex","first","cdid","mktsize","popular"]
            k_sym = (k == "first") ? :d_first : Symbol(k)
	        data[d_k][k_sym] = vecI64(dataset[k])
	    end
		# Type cast floats
		for k in ["numlist","localint","conditiondif","p","obsweight","condition"]
            k_sym = (k == "conditiondif") ? :cond_dif : (k == "obsweight") ? :obs_w : Symbol(k)
	        data[d_k][k_sym] = vecF64(dataset[k])
	   	end

		if d_k != :of_09
			data[d_k][:basecond] = vecF64(dataset["basecond"])
			data[d_k][:pdif]     = vecF64(dataset["pdif"])
		end

		if d_k == :on_12
			data[d_k][:disappear] = vecF64(dataset["disappear"])
		end
	end
	@save "$INPUT/data_to_run.jld2" data
	return data
end
