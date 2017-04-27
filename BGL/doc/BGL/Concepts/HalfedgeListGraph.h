/*!
\ingroup PkgBGLConcepts
\cgalConcept

The concept `HalfedgeListGraph` refines the concept `HalfedgeGraph`
and adds the requirements for traversal of all halfedges in the graph.

\cgalRefines `HalfedgeGraph`
\cgalHasModel `CGAL::Polyhedron_3`
\cgalHasModel `CGAL::Surface_mesh`

*/

class HalfedgeListGraph {};

/*! \relates HalfedgeListGraph
 * An iterator range over all halfedges.
 */

std::pair<boost::graph_traits<HalfedgeListGraph>::halfedge_iterator,
          boost::graph_traits<HalfedgeListGraph>::halfedge_iterator>
halfedges(const HalfedgeListGraph& g);


/*! \relates HalfedgeListGraph
  returns an upper bound of the number of halfedges of the graph.
  \attention `num_halfedges()` may return a number larger than `std::distance(halfedges(g).first, halfedges(g).second)`.
  This is the case for implementations only marking halfedges deleted in the halfedge container.
 */

boost::graph_traits<HalfedgeListGraph>::halfedge_size_type
num_halfedges(const FaceListGraph& g);

