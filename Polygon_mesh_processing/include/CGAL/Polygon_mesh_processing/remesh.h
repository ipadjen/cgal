// Copyright (c) 2015 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMESH_H
#define CGAL_POLYGON_MESH_PROCESSING_REMESH_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/Uniform_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/Adaptive_sizing_field.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#ifdef CGAL_PMP_REMESHING_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

/*! \todo document me or merge the doc with the original overload*/
template<typename PolygonMesh
       , typename FaceRange
       , typename SizingFunction
       , typename NamedParameters>
void isotropic_remeshing(const FaceRange& faces
                       , SizingFunction& sizing
                       , PolygonMesh& pmesh
                       , const NamedParameters& np);

/*!
* \ingroup PMP_meshing_grp
*
* @brief remeshes a triangulated region of a polygon mesh.
* This operation sequentially performs edge splits, edge collapses,
* edge flips, tangential relaxation and projection to the initial surface
* to generate a smooth mesh with a prescribed edge length.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         The descriptor types `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         and `boost::graph_traits<PolygonMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param pmesh a polygon mesh with triangulated surface patches to be remeshed
* @param faces the range of triangular faces defining one or several surface patches to be remeshed
* @param target_edge_length the edge length that is targeted in the remeshed patch.
*        If `0` is passed then only the edge-flip, tangential relaxation, and projection steps will be done.
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* @pre if constraints protection is activated, the constrained edges must
* not be longer than 4/3*`target_edge_length`.
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*     \cgalParamExtra{Exact constructions kernels are not supported by this function.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{face_index_map}
*     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
*                    as key type and `std::size_t` as value type}
*     \cgalParamDefault{an automatically indexed internal map}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of iterations for the sequence of atomic operations performed (listed in the above description)}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{edge_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no edge is constrained}
*     \cgalParamExtra{A constrained edge can be split or collapsed, but not flipped, nor its endpoints moved by smoothing.}
*     \cgalParamExtra{Sub-edges generated by splitting are set to be constrained.}
*     \cgalParamExtra{Patch boundary edges (i.e. incident to only one face in the range) are always considered as constrained edges.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of `pmesh`.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no vertex is constrained}
*     \cgalParamExtra{A constrained vertex cannot be modified during remeshing.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{protect_constraints}
*     \cgalParamDescription{If `true`, the edges set as constrained in `edge_is_constrained_map`
*                           (or by default the boundary edges) are neither split nor collapsed during remeshing.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*     \cgalParamExtra{Note that around constrained edges that have their length higher than
*                     twice `target_edge_length`, remeshing will fail to provide good quality results.
*                     It can even fail to terminate because of cascading vertex insertions.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{collapse_constraints}
*     \cgalParamDescription{If `true`, the edges set as constrained in `edge_is_constrained_map`
*                           (or by default the boundary edges) are collapsed during remeshing.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*     \cgalParamExtra{This value is ignored if `protect_constraints` is `true`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{face_patch_map}
*     \cgalParamDescription{a property map with the patch id's associated to the faces of `faces`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
*                    as key type and the desired property, model of `CopyConstructible` and `LessThanComparable` as value type.}
*     \cgalParamDefault{a default property map where each face is associated with the ID of
*                       the connected component it belongs to. Connected components are
*                       computed with respect to the constrained edges listed in the property map
*                       `edge_is_constrained_map`}
*     \cgalParamExtra{The map is updated during the remeshing process while new faces are created.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{do_split}
*     \cgalParamDescription{whether edges that are too long with respect to the given sizing are split}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{do_collapse}
*     \cgalParamDescription{whether edges that are too short with respect to the given sizing are collapsed}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{do_flip}
*     \cgalParamDescription{whether edge flips are performed to improve shape and valence}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{number_of_relaxation_steps}
*     \cgalParamDescription{the number of iterations of tangential relaxation that are performed
*                           at each iteration of the remeshing process}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{relax_constraints}
*     \cgalParamDescription{If `true`, the end vertices of the edges set as constrained
*                           in `edge_is_constrained_map` and boundary edges move along the
*                           constrained polylines they belong to.}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`false`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{do_project}
*     \cgalParamDescription{whether vertices should be reprojected on the input surface after creation or displacement}
*     \cgalParamType{Boolean}
*     \cgalParamDefault{`true`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{projection_functor}
*     \cgalParamDescription{A function object used to project input vertices (moved by the smoothing) and created vertices}
*     \cgalParamType{Unary functor that provides `%Point_3 operator()(vertex_descriptor)`, `%Point_3` being the value type
*                    of the vertex point map.}
*     \cgalParamDefault{If not provided, vertices are projected on the input surface mesh.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @sa `split_long_edges()`
*
*@todo Deal with exact constructions Kernel. The only thing that makes sense is to
*      guarantee that the output vertices are exactly on the input surface.
*      To do so, we can do every construction in `double`, and use an exact process for
*      projection. For each vertex, the `AABB_tree` would be used in an inexact manner
*      to find the triangle on which projection has to be done. Then, use
*      `CGAL::intersection(triangle, line)` in the exact constructions kernel to
*      get a point which is exactly on the surface.
*
*/
template<typename PolygonMesh
       , typename FaceRange
       , typename NamedParameters = parameters::Default_named_parameters>
void isotropic_remeshing(const FaceRange& faces
                       , const double target_edge_length
                       , PolygonMesh& pmesh
                       , const NamedParameters& np = parameters::default_values())
{
  typedef Uniform_sizing_field<PolygonMesh> Default_sizing;
  Default_sizing sizing(target_edge_length, pmesh);
  isotropic_remeshing<PolygonMesh, FaceRange, Default_sizing, NamedParameters>(
    faces,
    sizing,
    pmesh,
    np);
}

//todo ip: should I have the overload here?
/*
template<typename PolygonMesh
        , typename FaceRange
        , typename NamedParameters = parameters::Default_named_parameters>
void isotropic_remeshing(const FaceRange& faces
                       , const double& tol
                       , const std::pair<double, double>& edge_len_min_max
                       , PolygonMesh& pmesh
                       , const NamedParameters& np = parameters::default_values())
{
  typedef Adaptive_sizing_field<PolygonMesh> Adaptive_sizing;
  Adaptive_sizing sizing(edge_len_min_max, pmesh);
  isotropic_remeshing<PolygonMesh, FaceRange, Adaptive_sizing, NamedParameters>(
          faces,
          sizing,
          pmesh,
          np);
}
 */

template<typename PolygonMesh
       , typename FaceRange
       , typename SizingFunction
       , typename NamedParameters>
void isotropic_remeshing(const FaceRange& faces
                       , SizingFunction& sizing
                       , PolygonMesh& pmesh
                       , const NamedParameters& np)
{
  if (boost::begin(faces)==boost::end(faces))
    return;

  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;

  using parameters::get_parameter;
  using parameters::choose_parameter;

#ifdef CGAL_PMP_REMESHING_VERBOSE
  std::cout << std::endl;
  CGAL::Timer t;
  std::cout << "Remeshing parameters...";
  std::cout.flush();
  t.start();
#endif

  static const bool need_aabb_tree =
    parameters::is_default_parameter<NamedParameters, internal_np::projection_functor_t>::value;

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_property_map(vertex_point, pmesh));

  typedef typename GetInitializedFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
  FIMap fimap = CGAL::get_initialized_face_index_map(pmesh, np);

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::edge_is_constrained_t,
      NamedParameters,
      Static_boolean_property_map<edge_descriptor, false> // default (no constraint pmap)
    > ::type ECMap;
  ECMap ecmap = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                                 Static_boolean_property_map<edge_descriptor, false>());

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      Static_boolean_property_map<vertex_descriptor, false> // default (no constraint pmap)
    > ::type VCMap;
  VCMap vcmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                 Static_boolean_property_map<vertex_descriptor, false>());

  bool protect = choose_parameter(get_parameter(np, internal_np::protect_constraints), false);
  typedef typename internal_np::Lookup_named_param_def <
      internal_np::face_patch_t,
      NamedParameters,
      internal::Connected_components_pmap<PM, FIMap>//default
    > ::type FPMap;
  FPMap fpmap = choose_parameter(
    get_parameter(np, internal_np::face_patch),
    internal::Connected_components_pmap<PM, FIMap>(faces, pmesh, ecmap, fimap,
      parameters::is_default_parameter<NamedParameters, internal_np::face_patch_t>::value && (need_aabb_tree
#if !defined(CGAL_NO_PRECONDITIONS)
      || protect // face patch map is used to identify patch border edges to check protected edges are short enough
#endif
    ) ) );

#if !defined(CGAL_NO_PRECONDITIONS)
  if(protect)
  {
    std::string msg("Isotropic remeshing : protect_constraints cannot be set to");
    msg.append(" true with constraints larger than 4/3 * target_edge_length.");
    msg.append(" Remeshing aborted.");
    CGAL_precondition_msg(
      internal::constraints_are_short_enough(pmesh, ecmap, vpmap, fpmap, sizing),
      msg.c_str());
  }
#endif

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "\rRemeshing parameters done ("<< t.time() <<" sec)" << std::endl;
  std::cout << "Remesher construction...";
  std::cout.flush();
  t.reset(); t.start();
#endif

  typename internal::Incremental_remesher<PM, VPMap, GT, ECMap, VCMap, FPMap, FIMap>
    remesher(pmesh, vpmap, gt, protect, ecmap, vcmap, fpmap, fimap, need_aabb_tree);
  remesher.init_remeshing(faces);

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
#endif

  bool collapse_constraints = choose_parameter(get_parameter(np, internal_np::collapse_constraints), true);
  unsigned int nb_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);
  bool smoothing_1d = choose_parameter(get_parameter(np, internal_np::relax_constraints), false);
  unsigned int nb_laplacian = choose_parameter(get_parameter(np, internal_np::number_of_relaxation_steps), 1);
  bool do_collapse = choose_parameter(get_parameter(np, internal_np::do_collapse), true);
  bool do_split = choose_parameter(get_parameter(np, internal_np::do_split), true);
  bool do_flip = choose_parameter(get_parameter(np, internal_np::do_flip), true);

#ifdef CGAL_PMP_REMESHING_VERBOSE
  std::cout << std::endl;
  std::cout << "Remeshing (#iter = " << nb_iterations << ")..." << std::endl;
  t.reset(); t.start();
#endif

//      sizing.calc_sizing_map();
  for (unsigned int i = 0; i < nb_iterations; ++i)
  {
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

    if (i < 2)
      sizing.calc_sizing_map();
    if(do_split)
     remesher.split_long_edges(sizing);
    if(do_collapse)
     remesher.collapse_short_edges(sizing, collapse_constraints);
    if(do_flip)
      remesher.flip_edges_for_valence_and_shape();
    remesher.tangential_relaxation_impl(smoothing_1d, nb_laplacian);
    if ( choose_parameter(get_parameter(np, internal_np::do_project), true) )
      remesher.project_to_surface(get_parameter(np, internal_np::projection_functor));
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << std::endl;
#endif
  }

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "Remeshing done (#iter = " << nb_iterations;
  std::cout << ", " << t.time() << " sec )." << std::endl;
#endif
}

/*!
* \ingroup PMP_meshing_grp
* @brief splits the edges listed in `edges` into sub-edges
* that are not longer than the given threshold `max_length`.
*
* Note this function is useful to split constrained edges before
* calling `isotropic_remeshing()` with protection of constraints
* activated (to match the constrained edge length required by the
* remeshing algorithm to be guaranteed to terminate)
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`.
* @tparam EdgeRange range of `boost::graph_traits<PolygonMesh>::%edge_descriptor`,
*   model of `Range`. Its iterator type is `InputIterator`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param pmesh a polygon mesh
* @param edges the range of edges to be split if they are longer than given threshold
* @param max_length the edge length above which an edge from `edges` is split
*        into to sub-edges
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `PolygonMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{face_index_map}
*     \cgalParamDescription{a property map associating to each face of `pmesh` a unique index between `0` and `num_faces(pmesh) - 1`}
*     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
*                    as key type and `std::size_t` as value type}
*     \cgalParamDefault{an automatically indexed internal map}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{face_patch_map}
*     \cgalParamDescription{a property map with the patch id's associated to the faces of `faces`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor`
*                    as key type and the desired property, model of `CopyConstructible` and `LessThanComparable` as value type.}
*     \cgalParamDefault{a default property map where each face is associated with the ID of
*                       the connected component it belongs to. Connected components are
*                       computed with respect to the constrained edges listed in the property map
*                       `edge_is_constrained_map`}
*     \cgalParamExtra{The map is updated during the remeshing process while new faces are created.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{edge_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `pmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no edge is constrained}
*     \cgalParamExtra{A constrained edge can be split or collapsed, but not flipped, nor its endpoints moved by smoothing.}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* @sa `isotropic_remeshing()`
*
*/
template<typename PolygonMesh
       , typename EdgeRange
       , typename NamedParameters = parameters::Default_named_parameters>
void split_long_edges(const EdgeRange& edges
                    , const double& max_length
                    , PolygonMesh& pmesh
                    , const NamedParameters& np = parameters::default_values())
{
  typedef PolygonMesh PM;
  typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetGeomTraits<PM, NamedParameters>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename GetVertexPointMap<PM, NamedParameters>::type VPMap;
  VPMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                 get_property_map(vertex_point, pmesh));

  typedef typename GetInitializedFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
  FIMap fimap = CGAL::get_initialized_face_index_map(pmesh, np);

  typedef typename internal_np::Lookup_named_param_def <
        internal_np::edge_is_constrained_t,
        NamedParameters,
        Static_boolean_property_map<edge_descriptor, false> // default (no constraint pmap)
      > ::type ECMap;
  ECMap ecmap = choose_parameter(get_parameter(np, internal_np::edge_is_constrained),
                                 Static_boolean_property_map<edge_descriptor, false>());

  typedef typename internal_np::Lookup_named_param_def <
      internal_np::face_patch_t,
      NamedParameters,
      internal::Connected_components_pmap<PM, FIMap>//default
    > ::type FPMap;
  FPMap fpmap = choose_parameter(
    get_parameter(np, internal_np::face_patch),
    internal::Connected_components_pmap<PM, FIMap>(faces(pmesh), pmesh, ecmap, fimap, false));

  typename internal::Incremental_remesher<PM, VPMap, GT, ECMap,
    Static_boolean_property_map<vertex_descriptor, false>, // no constraint pmap
    FPMap,FIMap
  >
    remesher(pmesh, vpmap, gt, false/*protect constraints*/, ecmap,
             Static_boolean_property_map<vertex_descriptor, false>(),
             fpmap,
             fimap,
             false/*need aabb_tree*/);

  remesher.split_long_edges(edges, max_length);
}

} //end namespace Polygon_mesh_processing
} //end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
