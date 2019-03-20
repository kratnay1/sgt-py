import group_modules as gm

group_pairs = gm.get_space_subgroups(17,4)
pair1 = group_pairs[0]

supergroup = pair1.supergroup
subgroup = pair1.subgroup

complements = gm.get_sym_complements(17,4,supergroup.matrix, subgroup.size())
complement = complements[0]

product = subgroup.lin_rep * complement
