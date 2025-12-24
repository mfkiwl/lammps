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

#ifndef LMP_IMAGE_OBJECTS_H
#define LMP_IMAGE_OBJECTS_H

#include <array>
#include <vector>

namespace LAMMPS_NS {
class Image;
class Region;
namespace ImageObjects {
  constexpr int RESOLUTION = 36;    // default resolution for cylindrical objects
  constexpr int DEFLEVEL = 3;       // default refinement level for ellipsoids

  // custom data types for positions and triangles
  using vec3 = std::array<double, 3>;
  using triangle = std::array<vec3, 3>;

  class ArrowObj {
   public:
    // build an arrow template with length 1 in (1.0, 0.0, 0.0) direction as list of triangles
    void construct(double _tipl = 0.2, double _tipw = 0.1, double radius = 0.1,
                   int res = RESOLUTION);

    // construct a custom arrow from the template by scaling, rotating, and translating it.
    std::vector<triangle> transform(const vec3 &, const vec3 &, double, double);

    // draw custom arrow from unit template
    void draw(Image *, const double *, const double *, double, const double *, double, double);

   private:
    double tiplength;
    double tipwidth;
    double diameter;
    std::vector<triangle> triangles;
    int resolution;
  };

  class EllipsoidObj {
   public:
    // construct triangle mesh by refinining the triangles of an octahedron
    void construct(int level = DEFLEVEL);

    // draw ellipsoid from triangle mesh for ellipsoid particles
    void draw(Image *, int, const double *, const double *, const double *, const double *, double,
              double opacity = 1.0);

    // draw ellipsoid from triangle mesh for regions
    void draw(Image *, int, const double *, const double *, const double *, Region *, double,
              double opacity = 1.0);

   private:
    std::vector<triangle> triangles;

    void refine();
  };
}    // namespace ImageObjects
}    // namespace LAMMPS_NS

#endif
