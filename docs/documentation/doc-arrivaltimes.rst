Arrival Times Measurements
==========================

To perform measurements using the arrival times method rather than the Gaussian ACF methods, use the ``arrivaltimes`` module.

.. code-block:: python
    :linenos:

    import arrivaltimes

The ``arrivaltimes`` module is largely independent from the rest of FRBGUI and can be used solely through scripting. The function ``arrivaltimes.measureburst`` is the main entry point into the module and performs most of the work of measurement. The module can display an interactive plot window to help with marking bursts, masking channels, and determining initial guesses through ``measureburst``'s ``show`` option.


``arrivaltimes`` module
-----------------------

.. automodule:: arrivaltimes
    :members:
    :undoc-members:
