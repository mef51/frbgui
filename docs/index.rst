.. This file is a bit of a hack to make the first TOC entry expandable.
.. https://github.com/readthedocs/sphinx_rtd_theme/issues/445

FRBGui Documentation
====================

Welcome to FRBGui's documentation. FRBGui is a graphical user interface for measuring spectro-temporal properties of `Fast Radio Bursts`_ (FRBs) from their waterfalls using 2D autocorrelation functions (ACF). It can be used to split bursts with multiple components, change the dispersion measure (DM), add noise filters, and other preparation tasks before measurement. FRBGui is built with the `Dear PyGui <https://github.com/hoffstadt/DearPyGui>`_ interface library.

Click `here for a demo <https://www.youtube.com/watch?v=a5db2Ed7SFE>`_ of FRBGui.

Please visit :ref:`Overview and Installation <overview>` to get started.

.. _Fast Radio Bursts: https://en.wikipedia.org/wiki/Fast_radio_burst


.. toctree::
	:maxdepth: 2
	:caption: Introduction

	Overview & Installation <overview>

.. toctree::
	:maxdepth: 2
	:numbered:
	:caption: Tutorials

	tutorials/tut-preparing
	tutorials/tut-measuring

.. toctree::
	:maxdepth: 2
	:caption: Documentation

	documentation/doc-scripting

.. toctree::
	:maxdepth: 2
	:caption: Reference

	reference/ref-output
	reference/ref-acfmethod

