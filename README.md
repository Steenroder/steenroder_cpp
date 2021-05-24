# steenroder
Copyright © 2020–2021 Guillaume Tauzin.

### Description
``steenroder`` is a C++ library for the computation of persistent steenrod barcodes inspired by ``PHAT`` available under the GPLv3 license. It is part of a collaboration with Anibal Medina-Mardones and Umberto Lupo. See https://github.com/ammedmar/steenroder and https://github.com/ulupo/steenroder.

### Building

Flagser requires a C++14 compiler and cmake. Here is how to obtain, build, and
run flagser:

```sh
git clone git@gitlab.com:gtauzin/steenroder.git
cd steenroder
(mkdir -p build && cd build && cmake .. && make stn_double_2 && ./src/stn_double_2 ../examples/rp2.phat rp2_test)
```

### Contacts

gtauzin [at] protonmail.com
