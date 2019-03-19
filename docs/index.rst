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

And that's how you do that.


==============================

.. automodule:: group_modules
   :members: get_space_group, get_space_subgroups, loadGroup, loadCosetRep 

   .. autoclass:: LinRep 
       :members: 

   .. autoclass:: SpaceGroup
      :members: __init__, write_to_file
      
   .. autoclass:: SpaceGroupPair
      :members: __init__, is_normal



===============================

