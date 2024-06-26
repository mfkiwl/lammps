/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "nstencil_ghost_bin.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<int DIM_3D>
NStencilGhostBin<DIM_3D>::NStencilGhostBin(LAMMPS *lmp) : NStencil(lmp)
{
  xyzflag = 1;
}

/* ----------------------------------------------------------------------
   create stencil based on bin geometry and cutoff
------------------------------------------------------------------------- */

template<int DIM_3D>
void NStencilGhostBin<DIM_3D>::create()
{
  int i, j, k;

  nstencil = 0;

  for (k = -sz; k <= sz; k++) {
    for (j = -sy; j <= sy; j++) {
      for (i = -sx; i <= sx; i++) {
        if (bin_distance(i, j, k) < cutneighmaxsq) {
          stencilxyz[nstencil][0] = i;
          stencilxyz[nstencil][1] = j;
          stencilxyz[nstencil][2] = k;
          stencil[nstencil++] = k * mbiny * mbinx + j * mbinx + i;
        }
      }
    }
  }
}

namespace LAMMPS_NS {
template class NStencilGhostBin<0>;
template class NStencilGhostBin<1>;
}
