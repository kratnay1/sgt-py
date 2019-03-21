.. Sgt-Py documentation master file, created by
   sphinx-quickstart on Tue Nov 20 17:25:33 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Sgt-Py's documentation!

==================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

==================

* :ref:`genindex`
* :ref:`search`

===============================

Example usage:
===============================

Print the general positions of a space group::

   >>> import group_modules as gm
   >>> group_19 = gm.get_space_group(19)
   >>> print(group_19.lin_rep)
   x y z
   -x+1/2 -y z+1/2
   -x y+1/2 -z+1/2
   x+1/2 -y+1/2 -z

**Decomposing a space group into a product of subgroups**

To find a decomposition of the space group :math:`P222_1` with respect to the Bieberbach subgroup :math:`P2_1`, start by finding a SpaceGroupPair object relating the two space groups::

   >>> import group_modules as gm
   >>> group_pairs = gm.get_space_subgroups(17,4)
   >>> pair1 = group_pairs[0]
   >>> supergroup = pair1.supergroup
   >>> subgroup = pair1.subgroup

Find a complementary symmorphic subgroup of the right size::

   >>> complements = gm.get_sym_complements(17, 4, supergroup.matrix, subgroup.size())
   >>> complement = complements[0]

Check if the product of the two subgroups recovers the full group::

   >>> product = subgroup.lin_rep * complement
   >>> print(product == supergroup.lin_rep)
   True

Determine the type of product::

   >>> print(gm.is_normal(subgroup.lin_rep, supergroup.lin_rep))
   True
   >>> print(gm.is_normal(complement, supergroup.lin_rep))
   True

Identify the space group and transformation matrix that yields the complement group::

   >>> num, matrix = gm.identify_group(complement.cosets)
   >>> print(num)
   '3'

From this you can conclude that :math:`P222_1` can be written as a direct product of :math:`P2_1` and :math:`P3`.

==============================
Classes
==============================
.. automodule:: group_modules
   :members:

   .. autoclass:: LinRep 
       :members: 

   .. autoclass:: SpaceGroup
      :members: __init__, write_to_file
      
   .. autoclass:: SpaceGroupPair
      :members: __init__, is_normal


==============================
Functions
==============================

.. automodule:: group_modules
   :members: get_space_group, get_space_subgroups, loadGroup, loadCosetRep 



===============================

