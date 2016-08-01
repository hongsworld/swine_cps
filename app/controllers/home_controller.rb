class HomeController < ApplicationController
	def swine
	end
	def human
	end
	def chicken
	end

	def swine_submit
    @swine_aminoacid_full_name = @@aminoacid_full_name
    @swine_translation = @@translation
    swine_sample = params[:sequence].gsub(" ","")
    @swine_sequence_arr = transcription(swine_sample).scan(/.../)
    @swine_combination_cps = Hash.new
    swine_cpb_total = 0
    swine_cpb_under = Array.new
    @swine_cps_arr = Array.new
    @swine_minimal_candidate_arr = Array.new
    @swine_minimal_candidate_cpses = Array.new
    @swine_sequence_arr.each_with_index do |x,i|
      if i != (@swine_sequence_arr.count - 1) and x.index("N") == nil and @swine_sequence_arr[i+1].index("N") == nil
        swine_first = x
        swine_second = @swine_sequence_arr[i+1]
        if "#{@@translation[swine_first]}-#{@@translation[swine_second]}".index("X") == nil
          swine_cps = @@cps_swine_hash["#{reverse_transcription(swine_first)}-#{reverse_transcription(swine_second)}"]
          @swine_combination_cps[reverse_transcription("#{swine_first}-#{swine_second}")] = swine_cps
          @swine_cps_arr << swine_cps
          swine_cpb_total += swine_cps
          swine_cpb_under << "#{swine_first}-#{swine_second}"
          swine_first_codon_candidates = @@translation.select{|key,value| value == @@translation[swine_first]}.keys
          swine_second_codon_candidates = @@translation.select{|key,value| value == @@translation[swine_second]}.keys
          swine_codon_combination_candiates = Array.new
          swine_first_codon_candidates.each do |x|
            swine_second_codon_candidates.each do |y|
              swine_codon_combination_candiates << reverse_transcription("#{x}-#{y}")
            end
          end
          swine_codon_combination_cpses = Array.new
          swine_codon_combination_candiates.each do |x|
            swine_codon_combination_cpses << @@cps_swine_hash[x]
          end
          swine_minimal_candidate_cps = swine_codon_combination_cpses.min
          swine_minimal_candidate = swine_codon_combination_candiates[swine_codon_combination_cpses.index(swine_minimal_candidate_cps)]
          @swine_minimal_candidate_arr << swine_minimal_candidate
          @swine_minimal_candidate_cpses << swine_minimal_candidate_cps
        end
      end
    end
    @swine_cpb_score = swine_cpb_total/(swine_cpb_under.uniq.count)
    @swine_cps_arr << 0
    @swine_positive_cps_min = @swine_cps_arr.sort[@swine_cps_arr.sort.index(0) + 1]
    @swine_positive_cps_min = 0 if @swine_positive_cps_max == nil
    @swine_positive_cps_max = @swine_cps_arr.sort[-1]
    @swine_negative_cps_min = @swine_cps_arr.sort[0]
    @swine_negative_cps_max = @swine_cps_arr.sort[@swine_cps_arr.sort.index(0) - 1]
    @swine_negative_cps_max = 0 if @swine_positive_cps_max == nil
	end

	def human_submit
    @human_aminoacid_full_name = @@aminoacid_full_name
    @human_translation = @@translation
    human_sample = params[:sequence].gsub(" ","")
    @human_sequence_arr = transcription(human_sample).scan(/.../)
    @human_combination_cps = Hash.new
    human_cpb_total = 0
    human_cpb_under = Array.new
    @human_cps_arr = Array.new
    @human_minimal_candidate_arr = Array.new
    @human_minimal_candidate_cpses = Array.new
    @human_sequence_arr.each_with_index do |x,i|
      if i != (@human_sequence_arr.count - 1) and x.index("N") == nil and @human_sequence_arr[i+1].index("N") == nil
        human_first = x
        human_second = @human_sequence_arr[i+1]
        if "#{@@translation[human_first]}-#{@@translation[human_second]}".index("X") == nil
          human_cps = @@cps_human_hash["#{reverse_transcription(human_first)}-#{reverse_transcription(human_second)}"]
          @human_combination_cps[reverse_transcription("#{human_first}-#{human_second}")] = human_cps
          @human_cps_arr << human_cps
          human_cpb_total += human_cps
          human_cpb_under << "#{human_first}-#{human_second}"
          human_first_codon_candidates = @@translation.select{|key,value| value == @@translation[human_first]}.keys
          human_second_codon_candidates = @@translation.select{|key,value| value == @@translation[human_second]}.keys
          human_codon_combination_candiates = Array.new
          human_first_codon_candidates.each do |x|
            human_second_codon_candidates.each do |y|
              human_codon_combination_candiates << reverse_transcription("#{x}-#{y}")
            end
          end
          human_codon_combination_cpses = Array.new
          human_codon_combination_candiates.each do |x|
            human_codon_combination_cpses << @@cps_human_hash[x]
          end
          human_minimal_candidate_cps = human_codon_combination_cpses.min
          human_minimal_candidate = human_codon_combination_candiates[human_codon_combination_cpses.index(human_minimal_candidate_cps)]
          @human_minimal_candidate_arr << human_minimal_candidate
          @human_minimal_candidate_cpses << human_minimal_candidate_cps
        end
      end
    end
    @human_cpb_score = human_cpb_total/(human_cpb_under.uniq.count)
    @human_cps_arr << 0
    @human_positive_cps_min = @human_cps_arr.sort[@human_cps_arr.sort.index(0) + 1]
    @human_positive_cps_min = 0 if @human_positive_cps_max == nil
    @human_positive_cps_max = @human_cps_arr.sort[-1]
    @human_negative_cps_min = @human_cps_arr.sort[0]
    @human_negative_cps_max = @human_cps_arr.sort[@human_cps_arr.sort.index(0) - 1]
    @human_negative_cps_max = 0 if @human_positive_cps_max == nil
	end

	def chicken_submit
    @chicken_aminoacid_full_name = @@aminoacid_full_name
    @chicken_translation = @@translation
    chicken_sample = params[:sequence].gsub(" ","")
    @chicken_sequence_arr = transcription(chicken_sample).scan(/.../)
    @chicken_combination_cps = Hash.new
    chicken_cpb_total = 0
    chicken_cpb_under = Array.new
    @chicken_cps_arr = Array.new
    @chicken_minimal_candidate_arr = Array.new
    @chicken_minimal_candidate_cpses = Array.new
    @chicken_sequence_arr.each_with_index do |x,i|
      if i != (@chicken_sequence_arr.count - 1) and x.index("N") == nil and @chicken_sequence_arr[i+1].index("N") == nil
        chicken_first = x
        chicken_second = @chicken_sequence_arr[i+1]
        if "#{@@translation[chicken_first]}-#{@@translation[chicken_second]}".index("X") == nil
          chicken_cps = @@cps_chicken_hash["#{reverse_transcription(chicken_first)}-#{reverse_transcription(chicken_second)}"]
          @chicken_combination_cps[reverse_transcription("#{chicken_first}-#{chicken_second}")] = chicken_cps
          @chicken_cps_arr << chicken_cps
          chicken_cpb_total += chicken_cps
          chicken_cpb_under << "#{chicken_first}-#{chicken_second}"
          chicken_first_codon_candidates = @@translation.select{|key,value| value == @@translation[chicken_first]}.keys
          chicken_second_codon_candidates = @@translation.select{|key,value| value == @@translation[chicken_second]}.keys
          chicken_codon_combination_candiates = Array.new
          chicken_first_codon_candidates.each do |x|
            chicken_second_codon_candidates.each do |y|
              chicken_codon_combination_candiates << reverse_transcription("#{x}-#{y}")
            end
          end
          chicken_codon_combination_cpses = Array.new
          chicken_codon_combination_candiates.each do |x|
            chicken_codon_combination_cpses << @@cps_chicken_hash[x]
          end
          chicken_minimal_candidate_cps = chicken_codon_combination_cpses.min
          chicken_minimal_candidate = chicken_codon_combination_candiates[chicken_codon_combination_cpses.index(chicken_minimal_candidate_cps)]
          @chicken_minimal_candidate_arr << chicken_minimal_candidate
          @chicken_minimal_candidate_cpses << chicken_minimal_candidate_cps
        end
      end
    end
    @chicken_cpb_score = chicken_cpb_total/(chicken_cpb_under.uniq.count)
    @chicken_cps_arr << 0
    @chicken_positive_cps_min = @chicken_cps_arr.sort[@chicken_cps_arr.sort.index(0) + 1]
    @chicken_positive_cps_min = 0 if @chicken_positive_cps_max == nil
    @chicken_positive_cps_max = @chicken_cps_arr.sort[-1]
    @chicken_negative_cps_min = @chicken_cps_arr.sort[0]
    @chicken_negative_cps_max = @chicken_cps_arr.sort[@chicken_cps_arr.sort.index(0) - 1]
    @chicken_negative_cps_max = 0 if @chicken_positive_cps_max == nil
	end


  def united_submit
    @swine_aminoacid_full_name = @@aminoacid_full_name
    @swine_translation = @@translation
    swine_sample = params[:sequence].gsub(" ","")
    @swine_sequence_arr = transcription(swine_sample).scan(/.../)
    @swine_combination_cps = Hash.new
    swine_cpb_total = 0
    swine_cpb_under = Array.new
    @swine_cps_arr = Array.new
    @swine_minimal_candidate_arr = Array.new
    @swine_minimal_candidate_cpses = Array.new
    @swine_sequence_arr.each_with_index do |x,i|
      if i != (@swine_sequence_arr.count - 1) and x.index("N") == nil and @swine_sequence_arr[i+1].index("N") == nil
        swine_first = x
        swine_second = @swine_sequence_arr[i+1]
        if "#{@@translation[swine_first]}-#{@@translation[swine_second]}".index("X") == nil
          swine_cps = @@cps_swine_hash["#{reverse_transcription(swine_first)}-#{reverse_transcription(swine_second)}"]
          @swine_combination_cps[reverse_transcription("#{swine_first}-#{swine_second}")] = swine_cps
          @swine_cps_arr << swine_cps
          swine_cpb_total += swine_cps
          swine_cpb_under << "#{swine_first}-#{swine_second}"
          swine_first_codon_candidates = @@translation.select{|key,value| value == @@translation[swine_first]}.keys
          swine_second_codon_candidates = @@translation.select{|key,value| value == @@translation[swine_second]}.keys
          swine_codon_combination_candiates = Array.new
          swine_first_codon_candidates.each do |x|
            swine_second_codon_candidates.each do |y|
              swine_codon_combination_candiates << reverse_transcription("#{x}-#{y}")
            end
          end
          swine_codon_combination_cpses = Array.new
          swine_codon_combination_candiates.each do |x|
            swine_codon_combination_cpses << @@cps_swine_hash[x]
          end
          swine_minimal_candidate_cps = swine_codon_combination_cpses.min
          swine_minimal_candidate = swine_codon_combination_candiates[swine_codon_combination_cpses.index(swine_minimal_candidate_cps)]
          @swine_minimal_candidate_arr << swine_minimal_candidate
          @swine_minimal_candidate_cpses << swine_minimal_candidate_cps
        end
      end
    end
    @swine_cpb_score = swine_cpb_total/(swine_cpb_under.uniq.count)
    @swine_cps_arr << 0
    @swine_positive_cps_min = @swine_cps_arr.sort[@swine_cps_arr.sort.index(0) + 1]
    @swine_positive_cps_min = 0 if @swine_positive_cps_max == nil
    @swine_positive_cps_max = @swine_cps_arr.sort[-1]
    @swine_negative_cps_min = @swine_cps_arr.sort[0]
    @swine_negative_cps_max = @swine_cps_arr.sort[@swine_cps_arr.sort.index(0) - 1]
    @swine_negative_cps_max = 0 if @swine_positive_cps_max == nil

    @human_aminoacid_full_name = @@aminoacid_full_name
    @human_translation = @@translation
    human_sample = params[:sequence].gsub(" ","")
    @human_sequence_arr = transcription(human_sample).scan(/.../)
    @human_combination_cps = Hash.new
    human_cpb_total = 0
    human_cpb_under = Array.new
    @human_cps_arr = Array.new
    @human_minimal_candidate_arr = Array.new
    @human_minimal_candidate_cpses = Array.new
    @human_sequence_arr.each_with_index do |x,i|
      if i != (@human_sequence_arr.count - 1) and x.index("N") == nil and @human_sequence_arr[i+1].index("N") == nil
        human_first = x
        human_second = @human_sequence_arr[i+1]
        if "#{@@translation[human_first]}-#{@@translation[human_second]}".index("X") == nil
          human_cps = @@cps_human_hash["#{reverse_transcription(human_first)}-#{reverse_transcription(human_second)}"]
          @human_combination_cps[reverse_transcription("#{human_first}-#{human_second}")] = human_cps
          @human_cps_arr << human_cps
          human_cpb_total += human_cps
          human_cpb_under << "#{human_first}-#{human_second}"
          human_first_codon_candidates = @@translation.select{|key,value| value == @@translation[human_first]}.keys
          human_second_codon_candidates = @@translation.select{|key,value| value == @@translation[human_second]}.keys
          human_codon_combination_candiates = Array.new
          human_first_codon_candidates.each do |x|
            human_second_codon_candidates.each do |y|
              human_codon_combination_candiates << reverse_transcription("#{x}-#{y}")
            end
          end
          human_codon_combination_cpses = Array.new
          human_codon_combination_candiates.each do |x|
            human_codon_combination_cpses << @@cps_human_hash[x]
          end
          human_minimal_candidate_cps = human_codon_combination_cpses.min
          human_minimal_candidate = human_codon_combination_candiates[human_codon_combination_cpses.index(human_minimal_candidate_cps)]
          @human_minimal_candidate_arr << human_minimal_candidate
          @human_minimal_candidate_cpses << human_minimal_candidate_cps
        end
      end
    end
    @human_cpb_score = human_cpb_total/(human_cpb_under.uniq.count)
    @human_cps_arr << 0
    @human_positive_cps_min = @human_cps_arr.sort[@human_cps_arr.sort.index(0) + 1]
    @human_positive_cps_min = 0 if @human_positive_cps_max == nil
    @human_positive_cps_max = @human_cps_arr.sort[-1]
    @human_negative_cps_min = @human_cps_arr.sort[0]
    @human_negative_cps_max = @human_cps_arr.sort[@human_cps_arr.sort.index(0) - 1]
    @human_negative_cps_max = 0 if @human_positive_cps_max == nil

    @chicken_aminoacid_full_name = @@aminoacid_full_name
    @chicken_translation = @@translation
    chicken_sample = params[:sequence].gsub(" ","")
    @chicken_sequence_arr = transcription(chicken_sample).scan(/.../)
    @chicken_combination_cps = Hash.new
    chicken_cpb_total = 0
    chicken_cpb_under = Array.new
    @chicken_cps_arr = Array.new
    @chicken_minimal_candidate_arr = Array.new
    @chicken_minimal_candidate_cpses = Array.new
    @chicken_sequence_arr.each_with_index do |x,i|
      if i != (@chicken_sequence_arr.count - 1) and x.index("N") == nil and @chicken_sequence_arr[i+1].index("N") == nil
        chicken_first = x
        chicken_second = @chicken_sequence_arr[i+1]
        if "#{@@translation[chicken_first]}-#{@@translation[chicken_second]}".index("X") == nil
          chicken_cps = @@cps_chicken_hash["#{reverse_transcription(chicken_first)}-#{reverse_transcription(chicken_second)}"]
          @chicken_combination_cps[reverse_transcription("#{chicken_first}-#{chicken_second}")] = chicken_cps
          @chicken_cps_arr << chicken_cps
          chicken_cpb_total += chicken_cps
          chicken_cpb_under << "#{chicken_first}-#{chicken_second}"
          chicken_first_codon_candidates = @@translation.select{|key,value| value == @@translation[chicken_first]}.keys
          chicken_second_codon_candidates = @@translation.select{|key,value| value == @@translation[chicken_second]}.keys
          chicken_codon_combination_candiates = Array.new
          chicken_first_codon_candidates.each do |x|
            chicken_second_codon_candidates.each do |y|
              chicken_codon_combination_candiates << reverse_transcription("#{x}-#{y}")
            end
          end
          chicken_codon_combination_cpses = Array.new
          chicken_codon_combination_candiates.each do |x|
            chicken_codon_combination_cpses << @@cps_chicken_hash[x]
          end
          chicken_minimal_candidate_cps = chicken_codon_combination_cpses.min
          chicken_minimal_candidate = chicken_codon_combination_candiates[chicken_codon_combination_cpses.index(chicken_minimal_candidate_cps)]
          @chicken_minimal_candidate_arr << chicken_minimal_candidate
          @chicken_minimal_candidate_cpses << chicken_minimal_candidate_cps
        end
      end
    end
    @chicken_cpb_score = chicken_cpb_total/(chicken_cpb_under.uniq.count)
    @chicken_cps_arr << 0
    @chicken_positive_cps_min = @chicken_cps_arr.sort[@chicken_cps_arr.sort.index(0) + 1]
    @chicken_positive_cps_min = 0 if @chicken_positive_cps_max == nil
    @chicken_positive_cps_max = @chicken_cps_arr.sort[-1]
    @chicken_negative_cps_min = @chicken_cps_arr.sort[0]
    @chicken_negative_cps_max = @chicken_cps_arr.sort[@chicken_cps_arr.sort.index(0) - 1]
    @chicken_negative_cps_max = 0 if @chicken_positive_cps_max == nil
  end
end
