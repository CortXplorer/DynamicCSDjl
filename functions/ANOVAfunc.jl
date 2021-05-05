# function preANOVA(levels, variables)

# ID matrix should have in the first column the group identifier and in the second column the measurement identifier
id_matrix = Int.(ones(prod(levels),levels[1]))
id_matrix[:,1] = repeat(1:levels[1], inner=levels[2])
id_matrix[:,2] = repeat(1:levels[2], outer=levels[1])



# return X1, factorStarts, nColsFactor, labels, cellsets