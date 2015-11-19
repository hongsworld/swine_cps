triplet_codon_count = Hash.new
triplet_codon = Hash.new
triplet_codon["UUU"] = "F"
triplet_codon["UCU"] = "S"
triplet_codon["UAU"] = "Y"
triplet_codon["UGU"] = "C"
triplet_codon["UUC"] = "F"
triplet_codon["UCC"] = "S"
triplet_codon["UAC"] = "Y"
triplet_codon["UGC"] = "S"
triplet_codon["UUA"] = "L"
triplet_codon["UCA"] = "S"
triplet_codon["UAA"] = "X"
triplet_codon["UGA"] = "X"
triplet_codon["UUG"] = "L"
triplet_codon["UCG"] = "S"
triplet_codon["UAG"] = "X"
triplet_codon["UGG"] = "W"
triplet_codon["CUU"] = "L"
triplet_codon["CCU"] = "P"
triplet_codon["CAU"] = "H"
triplet_codon["CGU"] = "R"
triplet_codon["CUC"] = "L"
triplet_codon["CCC"] = "P"
triplet_codon["CAC"] = "H"
triplet_codon["CGC"] = "R"
triplet_codon["CUA"] = "L"
triplet_codon["CCA"] = "P"
triplet_codon["CAA"] = "Q"
triplet_codon["CGA"] = "R"
triplet_codon["CUG"] = "L"
triplet_codon["CCG"] = "P"
triplet_codon["CAG"] = "Q"
triplet_codon["CGG"] = "R"
triplet_codon["AUU"] = "I"
triplet_codon["ACU"] = "T"
triplet_codon["AAU"] = "N"
triplet_codon["AGU"] = "S"
triplet_codon["AUC"] = "I"
triplet_codon["ACC"] = "T"
triplet_codon["AAC"] = "N"
triplet_codon["AGC"] = "S"
triplet_codon["AUA"] = "I"
triplet_codon["ACA"] = "T"
triplet_codon["AAA"] = "K"
triplet_codon["AGA"] = "R"
triplet_codon["AUG"] = "M"
triplet_codon["ACG"] = "T"
triplet_codon["AAG"] = "K"
triplet_codon["AGG"] = "R"
triplet_codon["GUU"] = "V"
triplet_codon["GCU"] = "A"
triplet_codon["GAU"] = "D"
triplet_codon["GGU"] = "G"
triplet_codon["GUC"] = "V"
triplet_codon["GCC"] = "A"
triplet_codon["GAC"] = "D"
triplet_codon["GGC"] = "G"
triplet_codon["GUA"] = "V"
triplet_codon["GCA"] = "A"
triplet_codon["GAA"] = "E"
triplet_codon["GGA"] = "G"
triplet_codon["GUG"] = "V"
triplet_codon["GCG"] = "A"
triplet_codon["GAG"] = "E"
triplet_codon["GGG"] = "G"
#세팅하는것
triplet_codon.each do |a,b|
  triplet_codon_count[a] = 0
end

#세팅하는것
triplet_codon.each do |i,j|
  aminoacid_codon_pair[j] = Hash.new
end
#세팅하는것
triplet_codon.each do |i,j|
  aminoacid_codon_pair[j][i] = 0
end

codon_list = Array.new

triplet_codon.each do |a,b|
  codon_list << a
end

sample = %w[ATGTCGCCCATTCTCCTCTTCCATTCACCCCTTTCCTTTCTCTCCACAGCATTTCCAC
TGCTGGGCACGG
ATGTGGTCATGCACCCTGGGGCCTCACGCTCCTATCCTCCCCCCACTCCTGGCTCCCAACCCCAGCCAGC
TTGGGAGCAGCACCTGCCACCACCCATGGCCCCATTACTCCTGCCTGGCAGCACCCTATTGCAGCCAGCT
CACCCCATGACACCTCAGGTGACTCCTGCTGACCATGGCCCCATTGCCCCAGGGCTTGGCAACATCCCTG
TCCCAGTCTGGTCAGGAGGGAGGTCAGCATTGCCCCACCAGTCTCAGGCCTTTGGGCGCCCTCAGGCCCC
CCTCAATGGGAATGCTTCCAGTGCCCTCATGGGGAGTGCTTTGTGTCCTGCACCCCTTTGCATAGCAAGC
CCTGTGGGGGACAGCGGGGTGGCTGCCCCAGCTATTGCAGCAGCTCAGGCTGGCAAGGGAGGCTTGGCCC
CAGGCCTTCCTCCTCAAGCAGCACCGCCAGCTGCCCAGATGGCCTTCATCCTGCCCCCAGGGAGCTCTGG
GCCGTGGCCACATGGCACCTCTGGGGCAGGCAGCCTGCATGCCTCCCAGCACAAGGGCCCTCAGGAGGAC
GCCTCTCACCATAAGAGCGTCTATCAGAACTTCCGACGCTGGCAGCGCTTCAAGTCCCTGGCCCGCAGCC
ACCTGCCCCAGAGCCCCGATGCGGAAGCTCTCTCCTGCTTTCTCATCCCAGTGCTTCGCACCCTGGCTCG
CCTGAAGCCCACCATGACGCTGGAGGAGGGCATGCGGCAGGCCATGCAGGAATGGCAGCACATGAGCACC
GTTGAGAGGATGGCCTTTTATGAAATGGCAGCAAAGTTCATGCAATTTGAGGCAGAGGATGAGATCAGGA
TGGGGAACGTGTCCTGCACAACGTGGGCACAGGGCATGCCTCCTCCTGCTCCACCAAAGCTGGCTCCTCA
GGAACCCACATCCCGCAAAAGGGGCACTCAGTCTGTGAATCACCTTCCCCCTGGACGTCAAAGCACAGGG
GCCCCTGTGAAATCTGCCTAA]

def transcription(sample)
  return sample.gsub("A","1").gsub("T","2").gsub("C","3").gsub("G","4").gsub("1","U").gsub("2","A").gsub("3","G").gsub("4","C")
end



file_list = `ls swine_gene`.split("\n")

raw = Hash.new
each_sequence = Hash.new
total_sequence = Array.new
fail_sequence = Hash.new
file_list.each do |x|
  raw[x] = File.open("swine_gene/#{x}").read.gsub("\r","").split("\n");
  raw[x].each_with_index do |x,i|
    each_sequence[x] = x.split("\t")[1]
  end
  fail_sequence[x] = Array.new
  each_sequence[x].each do |y|
    if x.length%3 == 0
		total_sequence << y
	else
		fail_sequence[x] << y
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

