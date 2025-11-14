// clang-format off
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

#include "atom_vec_ellipsoid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm_kokkos.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory.h"
#include "memory_kokkos.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AtomVecEllipsoidKokkos::AtomVecEllipsoidKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecEllipsoid(lmp)
{
  //no_border_vel_flag = 0;
  //unpack_exchange_indices_flag = 1;
  size_border = 23;
  size_forward = 8;
  k_nghost_bonus = DAT::tdual_int_scalar("atomEllipKK:k_nghost_bonus");
  k_count_bonus = DAT::tdual_int_scalar("atomEllipKK:k_count_bonus");
}

/* ---------------------------------------------------------------------- */

AtomVecEllipsoidKokkos::~AtomVecEllipsoidKokkos()
{
  if (bonus_flag) {
    memoryKK->destroy_kokkos(k_bonus,bonus);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::init()
{
  AtomVecEllipsoid::init();

  set_atom_masks();
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by a chunk
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::grow(int n)
{
  auto DELTA = LMP_KOKKOS_AV_DELTA;
  int step = MAX(DELTA,nmax*0.01);
  if (n == 0) nmax += step;
  else nmax = n;
  atom->nmax = nmax;
  if (nmax < 0 || nmax > MAXSMALLINT)
    error->one(FLERR,"Per-processor system is too big");

  atomKK->sync(Device,ALL_MASK);
  atomKK->modified(Device,ALL_MASK);

  memoryKK->grow_kokkos(atomKK->k_tag,atomKK->tag,nmax,"atom:tag");
  memoryKK->grow_kokkos(atomKK->k_type,atomKK->type,nmax,"atom:type");
  memoryKK->grow_kokkos(atomKK->k_mask,atomKK->mask,nmax,"atom:mask");
  memoryKK->grow_kokkos(atomKK->k_image,atomKK->image,nmax,"atom:image");

  memoryKK->grow_kokkos(atomKK->k_x,atomKK->x,nmax,"atom:x");
  memoryKK->grow_kokkos(atomKK->k_v,atomKK->v,nmax,"atom:v");
  memoryKK->grow_kokkos(atomKK->k_f,atomKK->f,nmax,"atom:f");

  memoryKK->grow_kokkos(atomKK->k_rmass,atomKK->rmass,nmax,"atom:rmass");
  memoryKK->grow_kokkos(atomKK->k_angmom,atomKK->angmom,nmax,"atom:angmom");
  memoryKK->grow_kokkos(atomKK->k_torque,atomKK->torque,nmax,"atom:torque");
  memoryKK->grow_kokkos(atomKK->k_ellipsoid,atomKK->ellipsoid,nmax,"atom:ellipsoid");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);

  grow_pointers();
  atomKK->sync(Host,ALL_MASK);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::grow_pointers()
{
  tag = atomKK->tag;
  d_tag = atomKK->k_tag.view_device();
  h_tag = atomKK->k_tag.view_host();

  type = atomKK->type;
  d_type = atomKK->k_type.view_device();
  h_type = atomKK->k_type.view_host();
  mask = atomKK->mask;
  d_mask = atomKK->k_mask.view_device();
  h_mask = atomKK->k_mask.view_host();
  image = atomKK->image;
  d_image = atomKK->k_image.view_device();
  h_image = atomKK->k_image.view_host();

  x = atomKK->x;
  d_x = atomKK->k_x.view_device();
  h_x = atomKK->k_x.view_host();
  v = atomKK->v;
  d_v = atomKK->k_v.view_device();
  h_v = atomKK->k_v.view_host();
  f = atomKK->f;
  d_f = atomKK->k_f.view_device();
  h_f = atomKK->k_f.view_host();

  rmass = atomKK->rmass;
  d_rmass = atomKK->k_rmass.view_device();
  h_rmass = atomKK->k_rmass.view_host();
  angmom = atomKK->angmom;
  d_angmom = atomKK->k_angmom.view_device();
  h_angmom = atomKK->k_angmom.view_host();
  torque = atomKK->torque;
  d_torque = atomKK->k_torque.view_device();
  h_torque = atomKK->k_torque.view_host();
  ellipsoid = atomKK->ellipsoid;
  d_ellipsoid= atomKK->k_ellipsoid.view_device();
  h_ellipsoid = atomKK->k_ellipsoid.view_host();
}

/* ----------------------------------------------------------------------
   grow bonus data structure
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::grow_bonus()
{
  nmax_bonus = grow_nmax_bonus(nmax_bonus);
  if (nmax_bonus < 0) error->one(FLERR, "Per-processor system is too big");

  atomKK->sync(Device,BONUS_MASK);
  atomKK->modified(Device,BONUS_MASK);

  memoryKK->grow_kokkos(k_bonus,bonus,nmax_bonus,"atom:bonus");
  //d_bonus = k_bonus.d_view;
  //h_bonus = k_bonus.h_view;

  atomKK->sync(Host,BONUS_MASK);
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  atomKK->sync(Device, ALL_MASK & ~F_MASK & ~TORQUE_MASK);

  Sorter.sort(LMPDeviceType(), d_tag);
  Sorter.sort(LMPDeviceType(), d_type);
  Sorter.sort(LMPDeviceType(), d_mask);
  Sorter.sort(LMPDeviceType(), d_image);
  Sorter.sort(LMPDeviceType(), d_x);
  Sorter.sort(LMPDeviceType(), d_v);
  Sorter.sort(LMPDeviceType(), d_rmass);
  Sorter.sort(LMPDeviceType(), d_angmom);
  Sorter.sort(LMPDeviceType(), d_ellipsoid);
  Sorter.sort(LMPDeviceType(), d_bonus);

  atomKK->modified(Device, ALL_MASK & ~F_MASK & ~TORQUE_MASK);
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::sync(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.sync_device();
    if (mask & V_MASK) atomKK->k_v.sync_device();
    if (mask & F_MASK) atomKK->k_f.sync_device();
    if (mask & TAG_MASK) atomKK->k_tag.sync_device();
    if (mask & TYPE_MASK) atomKK->k_type.sync_device();
    if (mask & MASK_MASK) atomKK->k_mask.sync_device();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_device();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_device();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.sync_device();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_device();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.sync_device();
    if (mask & BONUS_MASK) k_bonus.sync_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.sync_host();
    if (mask & V_MASK) atomKK->k_v.sync_host();
    if (mask & F_MASK) atomKK->k_f.sync_host();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_host();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.sync_host();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_host();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.sync_host();
    if (mask & BONUS_MASK) k_bonus.sync_host();
  } else if (space == HostKK){
    if (mask & X_MASK) atomKK->k_x.sync_hostkk();
    if (mask & V_MASK) atomKK->k_v.sync_hostkk();
    if (mask & F_MASK) atomKK->k_f.sync_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.sync_host();
    if (mask & TYPE_MASK) atomKK->k_type.sync_host();
    if (mask & MASK_MASK) atomKK->k_mask.sync_host();
    if (mask & IMAGE_MASK) atomKK->k_image.sync_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.sync_hostkk();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.sync_hostkk();
    if (mask & TORQUE_MASK) atomKK->k_torque.sync_hostkk();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.sync_host();
    if (mask & BONUS_MASK) k_bonus.sync_host();
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag)
{
  if (space == Device) {
    if ((mask & X_MASK) && atomKK->k_x.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3_lr>(atomKK->k_x,space,async_flag);
    if ((mask & V_MASK) && atomKK->k_v.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_v,space,async_flag);
    if ((mask & F_MASK) && atomKK->k_f.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_f,space,async_flag);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync_device())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space,async_flag);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_type,space,async_flag);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_mask,space,async_flag);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync_device())
      perform_pinned_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space,async_flag);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & ANGMOM_MASK) && atomKK->k_angmom.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_angmom,space,async_flag);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync_device())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_torque,space,async_flag);
    if ((mask & ELLIPSOID_MASK) && atomKK->k_ellipsoid.need_sync_device())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_ellipsoid,space,async_flag);
    if ((mask & BONUS_MASK) && k_bonus.need_sync_device())
      perform_pinned_copy<DEllipsoidBonusAT::tdual_bonus_1d>(k_bonus,space,async_flag);
  } else {
    if ((mask & X_MASK) && atomKK->k_x.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3_lr>(atomKK->k_x,space,async_flag);
    if ((mask & V_MASK) && atomKK->k_v.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_v,space,async_flag);
    if ((mask & F_MASK) && atomKK->k_f.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_f,space,async_flag);
    if ((mask & TAG_MASK) && atomKK->k_tag.need_sync_host())
      perform_pinned_copy<DAT::tdual_tagint_1d>(atomKK->k_tag,space,async_flag);
    if ((mask & TYPE_MASK) && atomKK->k_type.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_type,space,async_flag);
    if ((mask & MASK_MASK) && atomKK->k_mask.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_mask,space,async_flag);
    if ((mask & IMAGE_MASK) && atomKK->k_image.need_sync_host())
      perform_pinned_copy<DAT::tdual_imageint_1d>(atomKK->k_image,space,async_flag);
    if ((mask & RMASS_MASK) && atomKK->k_rmass.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d>(atomKK->k_rmass,space,async_flag);
    if ((mask & ANGMOM_MASK) && atomKK->k_angmom.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_angmom,space,async_flag);
    if ((mask & TORQUE_MASK) && atomKK->k_torque.need_sync_host())
      perform_pinned_copy_transform<DAT::ttransform_kkfloat_1d_3>(atomKK->k_torque,space,async_flag);
    if ((mask & ELLIPSOID_MASK) && atomKK->k_ellipsoid.need_sync_host())
      perform_pinned_copy<DAT::tdual_int_1d>(atomKK->k_ellipsoid,space,async_flag);
    if ((mask & BONUS_MASK) && k_bonus.need_sync_host())
      perform_pinned_copy<DEllipsoidBonusAT::tdual_bonus_1d>(k_bonus,space,async_flag);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::modified(ExecutionSpace space, uint64_t mask)
{
  if (space == Device) {
    if (mask & X_MASK) atomKK->k_x.modify_device();
    if (mask & V_MASK) atomKK->k_v.modify_device();
    if (mask & F_MASK) atomKK->k_f.modify_device();
    if (mask & TAG_MASK) atomKK->k_tag.modify_device();
    if (mask & TYPE_MASK) atomKK->k_type.modify_device();
    if (mask & MASK_MASK) atomKK->k_mask.modify_device();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_device();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_device();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.modify_device();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_device();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.modify_device();
    if (mask & BONUS_MASK) k_bonus.modify_device();
  } else if (space == Host) {
    if (mask & X_MASK) atomKK->k_x.modify_host();
    if (mask & V_MASK) atomKK->k_v.modify_host();
    if (mask & F_MASK) atomKK->k_f.modify_host();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_host();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.modify_host();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_host();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.modify_host();
    if (mask & BONUS_MASK) k_bonus.modify_host();
  } else if (space == HostKK) {
    if (mask & X_MASK) atomKK->k_x.modify_hostkk();
    if (mask & V_MASK) atomKK->k_v.modify_hostkk();
    if (mask & F_MASK) atomKK->k_f.modify_hostkk();
    if (mask & TAG_MASK) atomKK->k_tag.modify_host();
    if (mask & TYPE_MASK) atomKK->k_type.modify_host();
    if (mask & MASK_MASK) atomKK->k_mask.modify_host();
    if (mask & IMAGE_MASK) atomKK->k_image.modify_host();
    if (mask & RMASS_MASK) atomKK->k_rmass.modify_hostkk();
    if (mask & ANGMOM_MASK) atomKK->k_angmom.modify_hostkk();
    if (mask & TORQUE_MASK) atomKK->k_torque.modify_hostkk();
    if (mask & ELLIPSOID_MASK) atomKK->k_ellipsoid.modify_host();
    if (mask & BONUS_MASK) k_bonus.modify_host();
  }
}

/* ------------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackCommBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_um _buf;
  typename AT::t_int_1d_const _list;
  typename AT::t_int_1d _ellipsoid;
  typename AtomVecEllipsoidKokkosBonusArray
          <DeviceType>::t_bonus_1d _bonus;

  AtomVecEllipsoidKokkos_PackCommBonus(
    const AtomKokkos* atomKK,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const typename DAT::tdual_double_2d_lr &buf,
    const typename DAT::tdual_int_1d &list):
    _buf(buf.view<DeviceType>()),
    _list(list.view<DeviceType>()),
    _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
    _bonus(bonus.view<DeviceType>()) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    if (_ellipsoid(j) >= 0) {
      /*_buf(i,4) = _bonus(_ellipsoid(j)).shape[0];
      _buf(i,5) = _bonus(_ellipsoid(j)).shape[1];
      _buf(i,6) = _bonus(_ellipsoid(j)).shape[2];*/
      _buf(i,4) = _bonus(_ellipsoid(j)).quat[0];
      _buf(i,5) = _bonus(_ellipsoid(j)).quat[1];
      _buf(i,6) = _bonus(_ellipsoid(j)).quat[2];
      _buf(i,7) = _bonus(_ellipsoid(j)).quat[3];
    }
  }
};

/* ------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::pack_comm_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                                          const DAT::tdual_double_2d_lr &buf)
{
  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->sync(HostKK,ELLIPSOID_MASK|BONUS_MASK);
    struct AtomVecEllipsoidKokkos_PackCommBonus<LMPHostType> f(atomKK,k_bonus,buf,list);
    Kokkos::parallel_for(n,f); 
  } else {
    atomKK->sync(Device,ELLIPSOID_MASK|BONUS_MASK);
    struct AtomVecEllipsoidKokkos_PackCommBonus<LMPDeviceType> f(atomKK,k_bonus,buf,list);
    Kokkos::parallel_for(n,f);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_UnpackCommBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_um _buf;
  int _first;
  typename AtomVecEllipsoidKokkosBonusArray
          <DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;

  AtomVecEllipsoidKokkos_UnpackCommBonus(
    const AtomKokkos* atomKK,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const typename DAT::tdual_double_2d_lr &buf,
    const int& first):
    _buf(buf.view<DeviceType>()),
    _first(first),
    _bonus(bonus.view<DeviceType>()),
    _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    if (_ellipsoid(i+_first) >= 0) {
      /*_bonus(_ellipsoid(i+_first)).shape[0] = _buf(i,4);
      _bonus(_ellipsoid(i+_first)).shape[1] = _buf(i,5);
      _bonus(_ellipsoid(i+_first)).shape[2] = _buf(i,6);*/
      _bonus(_ellipsoid(i+_first)).quat[0] = _buf(i,4);
      _bonus(_ellipsoid(i+_first)).quat[1] = _buf(i,5);
      _bonus(_ellipsoid(i+_first)).quat[2] = _buf(i,6);
      _bonus(_ellipsoid(i+_first)).quat[3] = _buf(i,7);
    }
  }
};

/* ---------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::unpack_comm_bonus_kokkos(
  const int &n, const int &first,
  const DAT::tdual_double_2d_lr &buf) {

  if (lmp->kokkos->forward_comm_on_host) {
    atomKK->modified(HostKK,ELLIPSOID_MASK|BONUS_MASK);
    struct AtomVecEllipsoidKokkos_UnpackCommBonus<LMPHostType> f(
      atomKK,k_bonus,buf,first);
    Kokkos::parallel_for(n,f);
  } else {
    atomKK->modified(Device,ELLIPSOID_MASK|BONUS_MASK);
    struct AtomVecEllipsoidKokkos_UnpackCommBonus<LMPDeviceType> f(
      atomKK,k_bonus,buf,first);
    Kokkos::parallel_for(n,f); 
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_PackBorderBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr _buf;
  const typename AT::t_int_1d_const _list;
  const typename AtomVecEllipsoidKokkosBonusArray
           <DeviceType>::t_bonus_1d _bonus;
  const typename AT::t_int_1d _ellipsoid;

  AtomVecEllipsoidKokkos_PackBorderBonus(
    const AtomKokkos* atomKK,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const typename AT::t_double_2d_lr &buf,
    const typename AT::t_int_1d_const &list):
    _buf(buf),
    _list(list),
    _bonus(bonus.view<DeviceType>()),
    _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    const int j = _list(i);
    const int j_bonus = _ellipsoid(j);
    if (j_bonus < 0) {
      _buf(i,7) = d_ubuf(j_bonus).d;
    } else {
      _buf(i,7) = d_ubuf(j_bonus).d;
      _buf(i,8) = _bonus(j_bonus).shape[0];
      _buf(i,9) = _bonus(j_bonus).shape[1];
      _buf(i,10) = _bonus(j_bonus).shape[2];
      _buf(i,11) = _bonus(j_bonus).quat[0];
      _buf(i,12) = _bonus(j_bonus).quat[1];
      _buf(i,13) = _bonus(j_bonus).quat[2];
      _buf(i,14) = _bonus(j_bonus).quat[3];
    }
  }
};

/* ---------------------------------------------------------------------- */
void AtomVecEllipsoidKokkos::pack_border_bonus_kokkos(
  int n, DAT::tdual_int_1d k_sendlist, 
  DAT::tdual_double_2d_lr buf, ExecutionSpace space)
{
  if (space==HostKK) {
    AtomVecEllipsoidKokkos_PackBorderBonus<LMPHostType> f(
      atomKK, k_bonus, buf.view_host(), k_sendlist.view_host());
    Kokkos::parallel_for(n,f);
  } else {
    AtomVecEllipsoidKokkos_PackBorderBonus<LMPDeviceType> f(
      atomKK, k_bonus, buf.view_device(), k_sendlist.view_device());
    Kokkos::parallel_for(n,f);
  }
}

/* ------------------------------------------------------------------------- */

template<class DeviceType>
struct AtomVecEllipsoidKokkos_UnpackBorderBonus {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  typename AT::t_double_2d_lr_const _buf;
  int _first;
  typename AtomVecEllipsoidKokkosBonusArray
           <DeviceType>::t_bonus_1d _bonus;
  typename AT::t_int_1d _ellipsoid;
  const int _nlocal_bonus;
  typename AT::t_int_scalar _nghost_bonus;

  AtomVecEllipsoidKokkos_UnpackBorderBonus(
    const AtomKokkos* atomKK,
    const typename DEllipsoidBonusAT::tdual_bonus_1d &bonus,
    const typename AT::t_double_2d_lr_const &buf,
    const int& first,
    const int &nlocal_bonus, 
    typename AT::tdual_int_scalar &nghost_bonus):
    _buf(buf),
    _first(first),
    _bonus(bonus.view<DeviceType>()),
    _ellipsoid(atomKK->k_ellipsoid.view<DeviceType>()),
    _nlocal_bonus(nlocal_bonus),
    _nghost_bonus(nghost_bonus.template view<DeviceType>()) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {

    int ellipID = static_cast<int> (d_ubuf(_buf(i,7)).i);

    if (ellipID == 0 ) {
      _ellipsoid(i+_first) = -1;
    } else {
      int j = _nlocal_bonus + Kokkos::atomic_fetch_add(&_nghost_bonus(),1);
      _bonus(j).shape[0] = _buf(i,8);
      _bonus(j).shape[1] = _buf(i,9); 
      _bonus(j).shape[2] = _buf(i,10);
      _bonus(j).quat[0] = _buf(i,11);
      _bonus(j).quat[1] = _buf(i,12);
      _bonus(j).quat[2] = _buf(i,13);
      _bonus(j).quat[3] = _buf(i,14);
      _bonus(j).ilocal = i+_first;
      _ellipsoid(i+_first) = j; 
    }
  }
};

/* ------------------------------------------------------------------------- */

void AtomVecEllipsoidKokkos::unpack_border_bonus_kokkos(const int &n, const int &first,
                            const DAT::tdual_double_2d_lr &buf,ExecutionSpace space) {
  while (first+n >= nmax) grow(0);
  while (n+nlocal_bonus+nghost_bonus >= nmax_bonus) grow_bonus();

  if (space==HostKK) {
    k_nghost_bonus.h_view() = nghost_bonus;
    struct AtomVecEllipsoidKokkos_UnpackBorderBonus<LMPHostType> f(
      atomKK, k_bonus, buf.view_host(), first,
      this->nlocal_bonus, k_nghost_bonus);
    Kokkos::parallel_for(n,f);
  } else {
    k_nghost_bonus.h_view() = nghost_bonus;
    k_nghost_bonus.modify<LMPHostType>();
    k_nghost_bonus.sync<LMPDeviceType>();
    struct AtomVecEllipsoidKokkos_UnpackBorderBonus<LMPDeviceType> f(
      atomKK, k_bonus, buf.view_device(), first,
      this->nlocal_bonus, k_nghost_bonus);
    Kokkos::parallel_for(n,f);
    k_nghost_bonus.modify<LMPDeviceType>();
    k_nghost_bonus.sync<LMPHostType>();
  }
  atomKK->modified(space,ELLIPSOID_MASK|BONUS_MASK);
}

// /* ---------------------------------------------------------------------- */

// int AtomVecEllipsoidKokkos::get_status_nlocal_bonus() {
//   return nlocal_bonus;
// }

// /* ---------------------------------------------------------------------- */

// void AtomVecEllipsoidKokkos::set_status_nlocal_bonus(int nlocal_bonus) {
//   this->nlocal_bonus = nlocal_bonus;
// }