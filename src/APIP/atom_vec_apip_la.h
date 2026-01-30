/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(apip/la,AtomVecApipLa);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_APIP_LA_H
#define LMP_ATOM_VEC_APIP_LA_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecApipLa : public AtomVec {
 public:
  AtomVecApipLa(class LAMMPS *);

  void grow_pointers() override;

  void force_clear(int, size_t) override;
  void data_atom_post(int) override;

 protected:
  double *apip_lambda, *apip_la_inp, *apip_la_avg, *apip_la_norm, *apip_e_fast, *apip_e_precise;
  int *apip_lambda_required;
};

#ifndef LMP_ATOM_APIP_LAMBDA_REQUIRED
#define LMP_ATOM_APIP_LAMBDA_REQUIRED
namespace ApipLambdaRequired {

  enum {
    UNKNOWN = 0,
    SIMPLE = 1 << 0,
    NO_SIMPLE = 1 << 1,
    COMPLEX = 1 << 2,
    NO_COMPLEX = 1 << 3,
  };
}    // namespace ApipLambdaRequired
#endif

}    // namespace LAMMPS_NS

#endif
#endif
