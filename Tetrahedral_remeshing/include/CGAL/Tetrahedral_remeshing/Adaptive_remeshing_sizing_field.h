// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
// File Description : Defines a sizing field adapted to a triangulation
//******************************************************************************

#ifndef CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H
#define CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Tetrahedral_remeshing/Sizing_field.h>


namespace CGAL
{
namespace Tetrahedral_remeshing
{
/**
 * @class Adaptive_remeshing_sizing_field
 */
template <typename Tr>
class Adaptive_remeshing_sizing_field
  : public Sizing_field<typename Tr::Geom_traits>
{
  // Types
  typedef typename Tr::Geom_traits              GT;
  typedef typename Tr::Geom_traits::Point_3     Bare_point;
  typedef typename Tr::Point                    Weighted_point;
  typedef typename GT::FT                       FT;

  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;

private:
  /// A cell that is used to accelerate location queries
  mutable Cell_handle last_cell_;

public:
  /**
   * Constructor
   */
  Adaptive_remeshing_sizing_field(Tr& tr);

  /**
   * Fills sizing field, using size associated to points in `tr_`
   */
  void fill();

  /**
  * Fills sizing field, using size associated to point in `value_map`
  */
  void fill(const std::map<Bare_point, FT>& value_map);

  /**
   * Returns size at point `p`.
   */
  FT operator()(const Bare_point& p) const
  { return this->operator()(p, this->get_last_cell()); }

  /**
   * Returns size at point `p`, using `v` to accelerate `p` location
   * in triangulation
   */
  FT operator()(const Bare_point& p, const Vertex_handle& v) const
  { return this->operator()(p, v->cell()); }

  /**
   * Returns size at point `p`.
   */
  FT operator()(const Bare_point& p, const Cell_handle& c) const;

  /**
   * Returns size at point `p`. Assumes that p is the centroid of c.
   */
  FT operator()(const Bare_point& p, const std::pair<Cell_handle, bool>& c) const;

  template <typename Index>
  FT operator()(const Bare_point& p, const int& dim, const Index& i) const
  { return this->operator()(p); }

  void set_triangulation(Tr& tr) { tr_ = tr; }

protected:
  Cell_handle get_last_cell() const
  {
    return last_cell_;
  }

  void set_last_cell(Cell_handle c) const
  {
    last_cell_ = c;
  }

private:
  /**
   * Returns size at point `p`, by interpolation into tetrahedron.
   */
  FT interpolate_on_cell_vertices(const Bare_point& p,
                                  const Cell_handle& cell) const;

  /**
   * Returns size at point `p`, by interpolation into facet (`cell` is assumed
   * to be an infinite cell).
   */
  FT interpolate_on_facet_vertices(const Bare_point& p,
                                   const Cell_handle& cell) const;

  FT sq_circumradius_length(const Cell_handle& cell, const Vertex_handle& v) const;
  FT average_circumradius_length(const Vertex_handle& v) const;

private:
  /// The triangulation
  Tr& tr_;
};


template <typename Tr>
Adaptive_remeshing_sizing_field<Tr>::
Adaptive_remeshing_sizing_field(Tr& tr)
  : tr_(tr)
{}

template <typename Tr>
void
Adaptive_remeshing_sizing_field<Tr>::
fill()
{
  std::map<Bare_point, FT> value_map;

  typename GT::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();

  // Fill map with local size
  for (typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
       vit != tr_.finite_vertices_end();
       ++vit)
  {
    const Weighted_point& position = tr_.point(vit);
    value_map.insert(std::make_pair(cp(position), average_circumradius_length(vit)));
  }

  // do the filling
  fill(value_map);
}

template <typename Tr>
void
Adaptive_remeshing_sizing_field<Tr>::
fill(const std::map<Bare_point, FT>& value_map)
{
  typedef typename Tr::Finite_vertices_iterator  Fvi;

  typename GT::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();

  for ( Fvi vit = tr_.finite_vertices_begin(); vit != tr_.finite_vertices_end(); ++ vit )
  {
    const Weighted_point& position = tr_.point(vit);

    typename std::map<Bare_point, FT>::const_iterator find_result =
      value_map.find(cp(position));

    if ( find_result != value_map.end() )
    {
      vit->set_meshing_info(find_result->second);
    }
    else
    {
      CGAL_assertion(false);
      vit->set_meshing_info(FT(0));
    }
  }
}


template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
operator()(const Bare_point& p, const Cell_handle& c) const
{
  typename GT::Construct_weighted_point_3 cwp = tr_.geom_traits().construct_weighted_point_3_object();

  //use the inexact locate (much faster than locate) to get a hint
  //and then use locate to check whether p is really inside hint
  // if not, an exact locate will be performed
  Cell_handle hint = tr_.inexact_locate(cwp(p),c);
  const Cell_handle cell = tr_.locate(cwp(p), hint);

  this->set_last_cell(cell);

  if ( !tr_.is_infinite(cell) )
    return interpolate_on_cell_vertices(p, cell);
  else
    return interpolate_on_facet_vertices(p, cell);
}


template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
operator()(const Bare_point&, const std::pair<Cell_handle,bool>& c) const
{
  // Assumes that p is the centroid of c
  const Cell_handle& cell = c.first;

  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(0)->meshing_info();
  const FT& vb = cell->vertex(1)->meshing_info();
  const FT& vc = cell->vertex(2)->meshing_info();
  const FT& vd = cell->vertex(3)->meshing_info();

  return ( (va+vb+vc+vd)/4 );
}


template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
interpolate_on_cell_vertices(const Bare_point& p, const Cell_handle& cell) const
{
  typename GT::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();
  typename GT::Compute_volume_3 volume = tr_.geom_traits().compute_volume_3_object();

  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(0)->meshing_info();
  const FT& vb = cell->vertex(1)->meshing_info();
  const FT& vc = cell->vertex(2)->meshing_info();
  const FT& vd = cell->vertex(3)->meshing_info();

  const Weighted_point& wa = tr_.point(cell, 0);
  const Weighted_point& wb = tr_.point(cell, 1);
  const Weighted_point& wc = tr_.point(cell, 2);
  const Weighted_point& wd = tr_.point(cell, 3);
  const Bare_point& a = cp(wa);
  const Bare_point& b = cp(wb);
  const Bare_point& c = cp(wc);
  const Bare_point& d = cp(wd);

  const FT abcp = CGAL::abs(volume(a,b,c,p));
  const FT abdp = CGAL::abs(volume(a,d,b,p));
  const FT acdp = CGAL::abs(volume(a,c,d,p));
  const FT bcdp = CGAL::abs(volume(b,d,c,p));

  // If volume is 0, then compute the average value
  if ( is_zero(abcp+abdp+acdp+bcdp) )
    return (va+vb+vc+vd) / 4.;

  return ( (abcp*vd + abdp*vc + acdp*vb + bcdp*va) / (abcp+abdp+acdp+bcdp) );
}


template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
interpolate_on_facet_vertices(const Bare_point& p, const Cell_handle& cell) const
{
  typename GT::Compute_area_3 area =  tr_.geom_traits().compute_area_3_object();

  typename GT::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();
  // Find infinite vertex and put it in k0
  int k0 = 0;
  int k1 = 1;
  int k2 = 2;
  int k3 = 3;

  if ( tr_.is_infinite(cell->vertex(1)) )
    std::swap(k0,k1);
  if ( tr_.is_infinite(cell->vertex(2)) )
    std::swap(k0,k2);
  if ( tr_.is_infinite(cell->vertex(3)) )
    std::swap(k0,k3);

  // Interpolate value using tet vertices values
  const FT& va = cell->vertex(k1)->meshing_info();
  const FT& vb = cell->vertex(k2)->meshing_info();
  const FT& vc = cell->vertex(k3)->meshing_info();

  const Weighted_point& wa = tr_.point(cell, k1);
  const Weighted_point& wb = tr_.point(cell, k2);
  const Weighted_point& wc = tr_.point(cell, k3);
  const Bare_point& a = cp(wa);
  const Bare_point& b = cp(wb);
  const Bare_point& c = cp(wc);

  const FT abp = area(a, b, p);
  const FT acp = area(a, c, p);
  const FT bcp = area(b, c, p);

  CGAL_assertion(abp >= 0);
  CGAL_assertion(acp >= 0);
  CGAL_assertion(bcp >= 0);

  // If area is 0, then compute the average value
  if ( is_zero(abp+acp+bcp) )
    return (va+vb+vc)/3.;

  return ( (abp*vc + acp*vb + bcp*va ) / (abp+acp+bcp) );
}

template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
sq_circumradius_length(const Cell_handle& cell, const Vertex_handle& v) const
{
  typename GT::Construct_point_3 cp
    = tr_.geom_traits().construct_point_3_object();
  typename GT::Compute_squared_distance_3 sq_distance
    = tr_.geom_traits().compute_squared_distance_3_object();
  typename GT::Construct_circumcenter_3 cc
    = tr_.geom_traits().construct_circumcenter_3_object();

  const auto t = tr_.tetrahedron(cell);
  const Bare_point circumcenter = cc(t[0], t[1], t[2], t[3]);
  const Weighted_point& position = tr_.point(cell, cell->index(v));

  return (sq_distance(cp(position), circumcenter));
}

template <typename Tr>
typename Adaptive_remeshing_sizing_field<Tr>::FT
Adaptive_remeshing_sizing_field<Tr>::
average_circumradius_length(const Vertex_handle& v) const
{
  std::vector<Cell_handle> incident_cells;
  incident_cells.reserve(64);
  tr_.incident_cells(v, std::back_inserter(incident_cells));

  using SI = typename Tr::Triangulation_data_structure::Cell::Subdomain_index;

  FT sum_len(0);
  unsigned int nb = 0;

  for (Cell_handle c : incident_cells)
  {
    if (c->subdomain_index() != SI())
    {
      sum_len += CGAL::approximate_sqrt(sq_circumradius_length(c, v));
      ++nb;
    }
  }

  // nb == 0 could happen if there is an isolated point.
  if (0 != nb)
  {
    return sum_len / nb;
  }
  else
  {
    // Use outside cells to compute size of point
    for (Cell_handle c : incident_cells)
    {
      if (!tr_.is_infinite(c))
      {
        sum_len += CGAL::approximate_sqrt(sq_circumradius_length(c, v));
        ++nb;
      }
    }

    CGAL_assertion(nb != 0);
    CGAL_assertion(sum_len != 0);
    return sum_len / nb;
  }
}

} // end namespace Tetrahedral_remeshing

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TETRAHEDRAL_REMESHING_ADAPTIVE_SIZING_FIELD_H
