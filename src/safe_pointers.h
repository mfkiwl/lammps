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

#ifndef LMP_SAFE_POINTERS_H
#define LMP_SAFE_POINTERS_H

// collection of smart pointers for specific purposes

#include <cstdio>

namespace LAMMPS_NS {

/** Class to automatically close a FILE pointer when a class instance goes out of scope.

\verbatim embed:rst

Drop in replacement for ``FILE *``. Use as ``SafeFilePtr fp;`` instead of
``FILE *fp = nullptr;`` and there is no more need to explicitly call
``fclose(fp)`` or ``pclose(fp)``.

\endverbatim
*/
class SafeFilePtr {
 public:
  SafeFilePtr() : fp(nullptr), use_pclose(false) {};
  SafeFilePtr(bool _use_pclose) : fp(nullptr), use_pclose(_use_pclose) {};
  SafeFilePtr(FILE *_fp, bool _use_pclose = false) : fp(_fp), use_pclose(_use_pclose) {};

  SafeFilePtr(const SafeFilePtr &) = delete;
  SafeFilePtr(SafeFilePtr &&o) noexcept : fp(o.fp), use_pclose(o.use_pclose) { o.fp = nullptr; }
  SafeFilePtr &operator=(const SafeFilePtr &) = delete;

  ~SafeFilePtr()
  {
    if (fp) {
      if (use_pclose)
        pclose(fp);
      else
        fclose(fp);
    }
  }

  SafeFilePtr &operator=(FILE *_fp)
  {
    if (fp && (fp != _fp)) {
      if (use_pclose)
        pclose(fp);
      else
        fclose(fp);
    }
    fp = _fp;
    return *this;
  }
  void set_pclose(bool _use_pclose) { use_pclose = _use_pclose; }
  operator FILE *() const { return fp; }

 private:
  FILE *fp;
  bool use_pclose;
};
}    // namespace LAMMPS_NS

#endif
