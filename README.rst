steenroder_cpp: A C++ library for the computation of persistence Steenrod barcodes
==================================================================================

Description
-----------
``steenroder_cpp`` is a C++ library for the computation of persistent steenrod barcodes inspired 
by ``PHAT`` available under the GPLv3 license. It was written by Guillaume Tauzin as part of a 
collaboration with Anibal Medina-Mardones and Umberto Lupo. See https://github.com/Steenroder/steenroder 
for a high-performance python package.

Context
-------
The widespread use in applied topology of the barcode of filtered
cellular complexes rests on a balance between discriminatory power and
computability. It has long been envision that the strength of this
invariant could be increase using cohomology operations. This package
computes the recently defined *Sq*\ \ *k*\ -barcodes which have been
shown to effectively increase the discriminatory power of barcodes on
real-world data.

For a complete presentation of these invariants please consult
`Persistence Steenrod modules <https://arxiv.org/abs/1812.05031>`__ by
U. Lupo, A. Medina-Mardones and G. Tauzin.

License
-------

``steenroder_cpp`` is distributed under the `GPLv3
license <https://github.com/Steenroder/steenroder_cpp/LICENSE>`__.

Building
--------

``steenroder_cpp`` requires a C++14 compiler and cmake. Here is how to obtain, build, and
run ``steenroder_cpp``:

```sh
git clone git@gitlab.com:Steenroder/steenroder_cpp.git
cd steenroder_cpp
(mkdir -p build && cd build && cmake .. && make stn_double_2 && ./src/stn_double_2 ../examples/rp2.phat rp2_test)
```

Important links
---------------

-  Official source code repo: https://github.com/Steenroder/steenroder_cpp
-  Issue tracker: https://github.com/Steenroder/steenroder_cpp/issues

Citing steenroder
-----------------

If you use ``steenroder_cpp`` in a scientific publication, we would
appreciate citations to the following paper:

`Persistence Steenrod modules <https://arxiv.org/abs/1812.05031>`__

You can use the following BibTeX entry:

::

   @article{steenroder,
          author = {{Lupo}, Umberto and {Medina-Mardones}, Anibal M. and {Tauzin}, Guillaume},
           title = "{Persistence Steenrod modules}",
         journal = {arXiv e-prints},
   archivePrefix = {arXiv},
          eprint = {1812.05031},
    primaryClass = {math.AT},
          adsurl = {https://ui.adsabs.harvard.edu/abs/2018arXiv181205031L},
   }

Contacts
--------

steenroder@gmail.com
