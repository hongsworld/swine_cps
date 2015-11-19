def transcription(sample)
	#return sample.gsub("A","1").gsub("T","2").gsub("C","3").gsub("G","4").gsub("1","U").gsub("2","A").gsub("3","G").gsub("4","C")
	return sample.gsub("T","U")
end
def reverse_transcription(sample)
	#return sample.gsub("A","1").gsub("T","2").gsub("C","3").gsub("G","4").gsub("1","U").gsub("2","A").gsub("3","G").gsub("4","C")
	return sample.gsub("U","T")
end
triplet_codon_count = Hash.new
translation = Hash.new
translation["UUU"] = "F"
translation["UCU"] = "S"
translation["UAU"] = "Y"
translation["UGU"] = "C"
translation["UUC"] = "F"
translation["UCC"] = "S"
translation["UAC"] = "Y"
translation["UGC"] = "S"
translation["UUA"] = "L"
translation["UCA"] = "S"
translation["UAA"] = "X"
translation["UGA"] = "X"
translation["UUG"] = "L"
translation["UCG"] = "S"
translation["UAG"] = "X"
translation["UGG"] = "W"
translation["CUU"] = "L"
translation["CCU"] = "P"
translation["CAU"] = "H"
translation["CGU"] = "R"
translation["CUC"] = "L"
translation["CCC"] = "P"
translation["CAC"] = "H"
translation["CGC"] = "R"
translation["CUA"] = "L"
translation["CCA"] = "P"
translation["CAA"] = "Q"
translation["CGA"] = "R"
translation["CUG"] = "L"
translation["CCG"] = "P"
translation["CAG"] = "Q"
translation["CGG"] = "R"
translation["AUU"] = "I"
translation["ACU"] = "T"
translation["AAU"] = "N"
translation["AGU"] = "S"
translation["AUC"] = "I"
translation["ACC"] = "T"
translation["AAC"] = "N"
translation["AGC"] = "S"
translation["AUA"] = "I"
translation["ACA"] = "T"
translation["AAA"] = "K"
translation["AGA"] = "R"
translation["AUG"] = "M"
translation["ACG"] = "T"
translation["AAG"] = "K"
translation["AGG"] = "R"
translation["GUU"] = "V"
translation["GCU"] = "A"
translation["GAU"] = "D"
translation["GGU"] = "G"
translation["GUC"] = "V"
translation["GCC"] = "A"
translation["GAC"] = "D"
translation["GGC"] = "G"
translation["GUA"] = "V"
translation["GCA"] = "A"
translation["GAA"] = "E"
translation["GGA"] = "G"
translation["GUG"] = "V"
translation["GCG"] = "A"
translation["GAG"] = "E"
translation["GGG"] = "G"
#세팅하는것
translation.each do |a,b|
	triplet_codon_count[a] = 0
end

#세팅하는것
aminoacid_codon_pair = Hash.new
translation.each do |a,b|
	aminoacid_codon_pair[b] = Hash.new
end
translation.each do |a,b|
	aminoacid_codon_pair[b][a] = 0
end
#세팅하는것
codon_list = Array.new
translation.each do |a,b|
	codon_list << a
end



file_list = `ls /Users/leehongseok/Developer/project/veterinary_computation/swine_project/swine_gene`.split("\n")

raw = Hash.new
each_sequence = Hash.new
each_sequence_info = Hash.new
total_gene_sequence_array = Array.new
fail_sequence = Hash.new
aminoacid_pool = Array.new
temp = 0
temp_success = 0
temp_3codon = 0

file_list.each do |x|
	puts x
  each_sequence[x] = Array.new
  each_sequence_info[x] = Array.new
  fail_sequence[x] = Array.new
  raw[x] = File.open("/Users/leehongseok/Developer/project/veterinary_computation/swine_project/swine_gene/#{x}").read.gsub("\r","").split("\n");
  raw[x].each_with_index do |z,i|
    each_sequence[x] << z.split("\t")[1]
	each_sequence_info[x] << z.split("\t")[0]
  end
  each_sequence[x].each_with_index do |y,i|
	temp += 1
	transcripted_sequence_codon = transcription(y).scan(/.../)
    if (y.length%3 == 0)
		temp_3codon += 1
		last_index = transcripted_sequence_codon.length - 1
		terminal_index_arr = transcripted_sequence_codon.map.with_index{|x, i| i if (x == "UGA") or (x == "UAA") or (x == "UAG")}.compact
		if  terminal_index_arr == [last_index] or terminal_index_arr == []
			total_gene_sequence_array <<  transcripted_sequence_codon
			temp_success += 1
		else
		  fail_sequence[x] << y
		end
    else
    end
  end
end

puts "총 sequence는 #{temp}"
puts "총 3codon sequence는 #{temp_3codon}"
puts "총 success sequence는 #{temp_success}"
#codon_combintation_count를 초기화합니다
codon_combination_count = Hash.new
codon_list.each do |x|
  codon_list.each do |y|
    codon_combination_count["#{x}-#{y}"] = 0
  end
end

aminoacid_combination_count = Hash.new
translation.each do |a,x|
	translation.each do |b,y|
		aminoacid_combination_count["#{x}-#{y}"] = 0
	end
end


#triplet_codon_count와 codon_combination_count, aminoacid_combination_count aminoacid_codon_pair_count를 얻습니다
total_gene_sequence_array.each do |each_gene_sequence|
	aminoacid_sequence = Array.new

	each_gene_sequence.each_with_index do |x,i|
		if x.index("N") == nil
			aminoacid_sequence << translation[x]
			triplet_codon_count[x] += 1

			if i != (each_gene_sequence.count - 1) and each_gene_sequence[i+1].index("N") == nil
				codon_combination_count["#{x}-#{each_gene_sequence[i+1]}"] += 1
				aminoacid_combination_count["#{translation[x]}-#{translation[each_gene_sequence[i+1]]}"] += 1
			end
		else
			#aminoacid_sequence에 더미를 하나 넣습니다
			aminoacid_sequence << ""
		end
	end

	aminoacid_sequence.each_with_index do |x,i|
		if x != ""
			if each_gene_sequence[i].index("N") == nil
				aminoacid_codon_pair[x][each_gene_sequence[i]] += 1
			end
		end
	end
end

#aminoacid_codon_pair bias를 계산합니다
aminoacid_codon_pair.each do |key,value|
	total = 0
	value.each do |a,b|
		total += b
	end
	if total != 0
		value.each do |a,b|
			value[a] = (b.to_f/total.to_f)
		end
	end
end

#combination_cps를 초기화합니다
combination_cps = Hash.new
codon_list.each do |x|
  codon_list.each do |y|
    combination_cps["#{x}-#{y}"] = 0
  end
end
puts aminoacid_combination_count
puts "\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
puts aminoacid_codon_pair
puts "\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

codon_combination_count.each do |a,count|
	x = a.split("-")[0]
	y = a.split("-")[1]
	if "#{translation[x]}-#{translation[y]}".index("X") == nil
		cps = Math.log( (count) / (aminoacid_combination_count["#{translation[x]}-#{translation[y]}"] * aminoacid_codon_pair[translation[x]][x] * aminoacid_codon_pair[translation[y]][y] ) )
		combination_cps[reverse_transcription(a)] = cps
	end
end


combination_cps.each do |x,y|
	puts "#{x}\t#{y}"
end
puts "\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
