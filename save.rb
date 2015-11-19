def transcription(sample)
	#return sample.gsub("A","1").gsub("T","2").gsub("C","3").gsub("G","4").gsub("1","U").gsub("2","A").gsub("3","G").gsub("4","C")
	return sample.gsub("T","U")
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
translation.each do |a,b|
	aminoacid_codon_pair[a] = Hash.new
	aminoacid_codon_pair[a][b] = 0
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

file_list.each do |x|
  each_sequence[x] = Array.new
  each_sequence_info[x] = Array.new
  fail_sequence[x] = Array.new
  raw[x] = File.open("/Users/leehongseok/Developer/project/veterinary_computation/swine_project/swine_gene/#{x}").read.gsub("\r","").split("\n");
  raw[x].each_with_index do |z,i|
    each_sequence[x] << z.split("\t")[1]
	each_sequence_info[x] << z.split("\t")[0]
  end
  each_sequence[x].each_with_index do |y,i|
	transcripted_sequence_codon = transcription(y).scan(/.../)
    if (y.length%3 == 0) and (transcripted_sequence_codon.index("UAA") != nil or transcripted_sequence_codon.index("UAG") != nil or transcripted_sequence_codon.index("UGA") != nil)
	  total_gene_sequence_array <<  transcripted_sequence_codon
    else
      fail_sequence[x] << y
    end
  end
end



sample = sample.join
#전사시킨다
sample = transcription(sample)

#세개씩 트리플렛 코돈으로 나눈다
codon = sample.scan(/.../)

amino_acid_sequence = Array.new

codon.each do |x|
	amino_acid_sequence << triplet_codon[x]
	triplet_codon_count[x] += 1
end

amino_acid_sequence.each_with_index do |x,i|
	aminoacid_codon_pair[x][codon[i]] += 1
end


#bias를 계산합니다
aminoacid_codon_pair.each do |key,value|
	total = 0
	value.each do |a,b|
		total += b
	end
	if total != 0
		value.each do |a,b|
			value[a] = (b.to_f/total.to_f)*(100.to_f)
		end
	end
end

puts aminoacid_codon_pair

codon_combination_count = Hash.new

codon_list.each do |x|
  codon_list.each do |y|
    codon_combination_count["#{x}-#{y}"] = 0
  end
end

codon.each_with_index do |x,i|
  if i != (codon.count - 1)
    codon_combination_count["#{x}-#{codon[i+1]}"] += 1
  end
end

combination_cps = Hash.new

codon_list.each do |x|
  codon_list.each do |y|
    combination_cps["#{triplet_codon[x]}-#{triplet_codon[y]}"] = 0
  end
end

codon_list.each do |x|
  codon_list.each do |y|
    cps = Math.log( (codon_combination_count["#{x}-#{y}"]) )
    combination_cps["#{triplet_codon[x]}-#{triplet_codon[y]}" = cps
  end
end

