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
* :ref:`modindex`
* :ref:`search`

===============================
===============================


.. code-block:: python
   import group_modules as gm

   grp19 = gm.get_spacegroup(19)
   print(grp19.lin_rep)

   grp_pairs = gm.get_space_subgroups(19,4)



.. automodule:: group_modules
   :members: get_space_group, get_space_subgroups, loadGroup, loadCosetRep 

   .. autoclass:: LinRep 
       :members: 

   .. autoclass:: SpaceGroup
      :members: __init__, write_to_file
      
   .. autoclass:: SpaceGroupPair
      :members: __init__, is_normal





===============================

