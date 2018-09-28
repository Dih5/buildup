=============================
buildup Python module
=============================

.. image:: https://readthedocs.org/projects/buildup/badge/?version=latest
        :target: https://buildup.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




A python module to use the build-up factors derived from methods described in the preprint available here_.

.. _here: https://arxiv.org/abs/1809.09907

.. image:: https://github.com/Dih5/buildup/raw/master/demos/demo.png


* Free software: LGPLv3+
* Documentation: https://buildup.readthedocs.io (available soon).


Features
--------
* Load existent build-up data for any response function (air exposure, tissue dose, ...).
* Obtain the uncertainty of the calculated factor due the numeric approximations.
* Generate and process FLUKA simulations to extend the database.

Instalation
-----------
To install the latest release, assuming you have a Python_ distribution with pip_::

    $ pip install https://github.com/Dih5/pileup/archive/master.zip
    
.. _Python: http://www.python.org/
.. _pip: https://pip.pypa.io/en/stable/installing/

This package includes the build-up description for the most commonly used materials. In the future, a method to download
only the desired materials from the applications will be added.


Usage
-----

Check the demo_ notebook.

.. _demo: https://github.com/Dih5/buildup/blob/master/demos/plotdemo.ipynb
