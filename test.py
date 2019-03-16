import group_modules as gm


subgroups = gm.get_space_subgroups(19,4,2)
p1 = subgroups[0]
sup1 = p1.supergroup
sub1 = p1.subgroup
g1 = sup1.lin_rep
s1 = sub1.lin_rep


p2 = subgroups[1]
sup2 = p2.supergroup
sub2 = p2.subgroup
g2 = sup2.lin_rep
s2 = sub2.lin_rep

