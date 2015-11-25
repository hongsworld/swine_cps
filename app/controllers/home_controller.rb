class HomeController < ApplicationController
	def index
	end

	def submit
		@aminoacid_full_name = @@aminoacid_full_name
		@translation = @@translation
		sample = params[:sequence].gsub(" ","")
		@sequence_arr = transcription(sample).scan(/.../)
		@combination_cps = Hash.new
		cpb_total = 0
		cpb_under = Array.new
		@cps_arr = Array.new
		@minimal_candidate_arr = Array.new
		@minimal_candidate_cpses = Array.new
		@sequence_arr.each_with_index do |x,i|
			if i != (@sequence_arr.count - 1) and x.index("N") == nil and @sequence_arr[i+1].index("N") == nil
				first = x
				second = @sequence_arr[i+1]
				if "#{@@translation[first]}-#{@@translation[second]}".index("X") == nil
					#cps = Math.log( (@sequence_arr.map.with_index{|x, i| i if (x == first) and (@sequence_arr[i+1] == second)}.compact.count) / (@@aminoacid_combination_count["#{@@translation[first]}-#{@@translation[second]}"] * @@aminoacid_codon_bias[@@translation[first]][first] * @@aminoacid_codon_bias[@@translation[second]][second] ) )
					cps = @@cps_hash["#{reverse_transcription(first)}-#{reverse_transcription(second)}"]
					@combination_cps[reverse_transcription("#{first}-#{second}")] = cps
					@cps_arr << cps
					cpb_total += cps
					cpb_under << "#{first}-#{second}"


					first_codon_candidates = @@translation.select{|key,value| value == @@translation[first]}.keys
					second_codon_candidates = @@translation.select{|key,value| value == @@translation[second]}.keys
					logger.info first_codon_candidates
					codon_combination_candiates = Array.new
					first_codon_candidates.each do |x|
						second_codon_candidates.each do |y|
							codon_combination_candiates << reverse_transcription("#{x}-#{y}")
						end
					end

					codon_combination_cpses = Array.new
					codon_combination_candiates.each do |x|
						codon_combination_cpses << @@cps_hash[x]
					end
					minimal_candidate_cps = codon_combination_cpses.min
					logger.info minimal_candidate_cps
					minimal_candidate = codon_combination_candiates[codon_combination_cpses.index(minimal_candidate_cps)]
					logger.info minimal_candidate
					@minimal_candidate_arr << minimal_candidate
					@minimal_candidate_cpses << minimal_candidate_cps
				end
			end
		end
		@cpb_score = cpb_total/(cpb_under.uniq.count)
		@cps_arr << 0
		@positive_cps_min = @cps_arr.sort[@cps_arr.sort.index(0) + 1]
		@positive_cps_min = 0 if @positive_cps_max == nil
		@positive_cps_max = @cps_arr.sort[-1]
		@negative_cps_min = @cps_arr.sort[0]
		@negative_cps_max = @cps_arr.sort[@cps_arr.sort.index(0) - 1]
		@negative_cps_max = 0 if @positive_cps_max == nil
	end
end
