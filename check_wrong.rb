def transcription(sample)
  return sample.gsub("T","U")
  #return sample.gsub("A","1").gsub("T","2").gsub("C","3").gsub("G","4").gsub("1","U").gsub("2","A").gsub("3","G").gsub("4","C").reverse
end

def reverse_transcription(sample)
  return sample.gsub("U","1").gsub("A","2").gsub("G","3").gsub("C","4").gsub("1","A").gsub("2","T").gsub("3","C").gsub("4","G").reverse
end

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
aminoacid_codon_pair = Hash.new
triplet_codon.each do |i,j|
  aminoacid_codon_pair[j] = Hash.new
end
#세팅하는것
triplet_codon.each do |i,j|
  aminoacid_codon_pair[j][i] = 0
end

triplet_codon_list = Array.new

triplet_codon.each do |a,b|
  triplet_codon_list << a
end

file_list = `ls /Users/leehongseok/Developer/project/veterinary_computation/swine_project/swine_gene`.split("\n")

raw = Hash.new
each_sequence = Hash.new
total_sequence = Array.new
fail_sequence = Hash.new
intact_gene_pool = Array.new
transcripted_gene_pool = Array.new
aminoacid_pool = Array.new
temp_check = 0
file_list.each do |x|
  raw[x] = File.open("/Users/leehongseok/Developer/project/veterinary_computation/swine_project/swine_gene/#{x}").read.gsub("\r","").split("\n");
  each_sequence[x] = Array.new
  raw[x].each_with_index do |z,i|
    each_sequence[x] << z.split("\t")[1]
  end
  fail_sequence[x] = Array.new
  each_sequence[x].each_with_index do |y,i|
    if y.length%3 == 0
      total_sequence << y
      if transcription(y).scan(/.../).index("UAA") != (y.scan(/.../).count - 1) and transcription(y).scan(/.../).index("UAA") != nil
		  temp_check += 1
		  puts "-------------------------------------"
		  puts transcription(y).scan(/.../).index("UAA")
		  puts "-------------------------------------"
		  puts "file"
		  puts x
		  puts "--------"
		  puts "index"
		  puts i
		  puts "--------"
		  puts "sequence"
		  puts y
	  end
    else
      fail_sequence[x] << y
    end
  end
end

puts temp_check
