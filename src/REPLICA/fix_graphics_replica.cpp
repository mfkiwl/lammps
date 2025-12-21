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

#include "fix_graphics_replica.h"

#include "atom.h"
#include "comm.h"
#include "dump_image.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::FixGraphicsReplica(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) error->universe_all(FLERR, "Too few arguments for fix graphics/replica");

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->universe_all(FLERR, "Illegal fix graphics/replica nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  dflag = false;
  aflag = false;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "display") == 0) {
      if (iarg + 4 > narg)
        error->universe_all(FLERR, "Too few arguments for fix graphics/replica display");
      dflag = true;
      iarg += 4;
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg + 4 > narg)
        error->universe_all(FLERR, "Too few arguments for fix graphics/replica average");
      aflag = true;
      iarg += 4;
    } else {
      error->universe_all(FLERR, std::string("Unknown fix graphics/replica keyword: ") + arg[iarg]);
    }
  }

  numobjs = 0;
}

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::~FixGraphicsReplica()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsReplica::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::init()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::end_of_step()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  bigint nper = 0;
  bigint nall = 0;
  if (comm->me == 0) nper = group->count(igroup);
  MPI_Allreduce(&nper, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, universe->uworld);

  bigint numtotal = 0;
  if (dflag) numtotal += nall;
  if (aflag) numtotal += nper;
  if (numtotal >= MAXSMALLINT) error->universe_all(FLERR, "Too many graphics objects");

  numobjs = (int) numtotal;

  if (universe->me == 0)
    fprintf(stderr, "We have %d atoms per replica and %d objects\n", nper, numobjs);

  memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
  memory->create(imgparms, numobjs, 5, "fix_graphics:imgparms");
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsReplica::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
