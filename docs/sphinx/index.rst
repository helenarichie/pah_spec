.. Here, we directly embed the contents of the README
   -> while it's okay right now to directly embed the entire README in this file
      we may change our mind in the future
   -> unforuntately, the myst_parser is fundamentally incompatible with including
      fragments that don't start with a `# header`. The parser is incompatible
      with starting with say a `## header`.
   -> so, if this ever comes up, we'll just need to duplicate contents OR
      convert the README to a .rst file (which GitHub definitely understands)

.. include:: ../../README.md
   :parser: myst_parser.sphinx_

.. Down below, we list contents (that populate the sidebar)

.. toctree::
   :maxdepth: 2
   :hidden:

   Introduction <self>
   install
   api
   development

.. toctree::
   :caption: Useful links
   :hidden:

   GitHub <https://github.com/helenarichie/pah_spec>
   Issue Tracker <https://github.com/helenarichie/pah_spec/issues>

.. Down below we describe information for the landing page of the website.

