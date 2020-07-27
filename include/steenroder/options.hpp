/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

namespace stn {

  enum class InputType { BoundaryMatrix, SimplexTree, DistanceMatrix, Image };
  enum class MatrixFormat { Generic, PHAT };
  enum class ComplexType { Simplicial, Cubical };
  enum class FieldType { Z2, Zp };

  enum class ReductionType { Generic, PHAT, Ripser, Flagser, GUDHI };


  enum class InputOutputFormat { Generic, None, PHAT, HDF5, XDMF };
  enum class InputOutputMode { Generic, ascii, binary };

}  // namespace stn
