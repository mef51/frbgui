Scripting
=========

FRBGui ships with the following modules that can be used for scripting outside of the graphical interface:

* :py:mod:`driftrate.py` : Utilities for manipulating FRB waterfalls and performing measurements on them.
* :py:mod:`driftlaw.py` : Utilities for processing FRBGui measurements for the purpose of studying spectro-temporal relationships.

These can be accessed with

.. code-block:: python
	:linenos:

	import driftrate
	import driftlaw

**This page is currently under construction and may be incomplete.**

``frbgui`` module
-----------------

the ``frbgui`` function itself has multiple options that can be used to configure the startup of FRBGui

As stated in the Overvew, FRBGui can be called from a script like so:

.. code-block:: python
	:linenos:

	from frbgui import frbgui

.. automodule:: frbgui
	:members:

``driftrate`` module
--------------------

This module contains Utilities for manipulating FRB waterfalls (such as dedispersion and subsampling) and performing measurements on them.

.. automodule:: driftrate
    :members:
    :undoc-members:

``driftlaw`` module
--------------------

.. automodule:: driftlaw
    :members:
    :undoc-members:
