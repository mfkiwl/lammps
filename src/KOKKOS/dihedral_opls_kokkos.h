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

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(opls/kk,DihedralOPLSKokkos<LMPDeviceType>);
DihedralStyle(opls/kk/device,DihedralOPLSKokkos<LMPDeviceType>);
DihedralStyle(opls/kk/host,DihedralOPLSKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_DIHEDRAL_OPLS_KOKKOS_H
#define LMP_DIHEDRAL_OPLS_KOKKOS_H

#include "dihedral_opls.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagDihedralOPLSCompute{};

template<class DeviceType>
class DihedralOPLSKokkos : public DihedralOPLS {
 public:
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;
  typedef ArrayTypes<DeviceType> AT;

  DihedralOPLSKokkos(class LAMMPS *);
  ~DihedralOPLSKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralOPLSCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagDihedralOPLSCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i1, const int i2, const int i3, const int i4,
                          F_FLOAT &edihedral, F_FLOAT *f1, F_FLOAT *f3, F_FLOAT *f4,
                          const F_FLOAT &vb1x, const F_FLOAT &vb1y, const F_FLOAT &vb1z,
                          const F_FLOAT &vb2x, const F_FLOAT &vb2y, const F_FLOAT &vb2z,
                          const F_FLOAT &vb3x, const F_FLOAT &vb3y, const F_FLOAT &vb3z) const;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;

 protected:

  class NeighborKokkos *neighborKK;
  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_2d dihedrallist;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_int_scalar k_warning_flag;
  typename AT::t_int_scalar d_warning_flag;
  HAT::t_int_scalar h_warning_flag;

  DAT::tdual_ffloat_1d k_k1;
  DAT::tdual_ffloat_1d k_k2;
  DAT::tdual_ffloat_1d k_k3;
  DAT::tdual_ffloat_1d k_k4;

  typename AT::t_ffloat_1d d_k1;
  typename AT::t_ffloat_1d d_k2;
  typename AT::t_ffloat_1d d_k3;
  typename AT::t_ffloat_1d d_k4;

  void allocate() override;
};

}

#endif
#endif

