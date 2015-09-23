=====================
Writing documentation
=====================

We use `Sphinx`_ for documentation.

We use the `autodoc`_ extension to generate source code documentation.

With the `Napoleon`_ extension, docstrings can be formatted nicely without the need to write RST formatting by hand.

With the `BibTeX`_ extension, we can use the same bibliography as in LaTeX.

.. _Sphinx: http://sphinx-doc.org/
.. _autodoc: http://sphinx-doc.org/ext/autodoc.html
.. _Napoleon: http://sphinx-doc.org/ext/napoleon.html
.. _BibTeX: http://sphinxcontrib-bibtex.readthedocs.org/


Compiling documentation
=======================

Needs:

* sphinx with some of the built-in extensions
* an external extension: sphinxcontrib-bibtex
* the Read the Docs theme

Almost certainly, these can all be installed from your favorite package
manager. If not, :code:`pip install --user ...` should always work.
