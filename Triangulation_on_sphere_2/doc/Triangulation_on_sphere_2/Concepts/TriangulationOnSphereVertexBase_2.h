/*!
\ingroup PkgTriangulationOnSphere2Concepts
\cgalConcept

\cgalRefines `TriangulationDSVertexBase_2`

The concept `TriangulationOnSphereVertexBase_2` describes the requirements for the vertex base class
of a triangulation data structure to be plugged in a Delaunay triangulation on the sphere.
It refines the concept `TriangulationDSVertexBase_2`, adding geometric information:
the vertex base of a triangulation stores a point.

\cgalHasModel `CGAL::Triangulation_on_sphere_vertex_base_2<Traits,Vb>`

\sa `TriangulationDataStructure_2`
\sa `TriangulationOnSphereFaceBase_2`
*/
class TriangulationOnSphereVertexBase_2
{
public:
  /// \name Types
  /// @{

  /*!
  Must be the equivalent to the point type `Point_on_sphere_2` of the concept `TriangulationOnSphereTraits_2`.
  */
  typedef unspecified_type Point;

  /*!
  Must be the equivalent to the point type `Face_handle` of the triangulation.
  */
  typedef unspecified_type Face_handle;

  /// @}

  /// \name Creation
  /// @{

  /*!
  Constructs a vertex whose geometric embedding is the point `p`.
  */
  TriangulationOnSphereVertexBase_2(Point p);

  /*!
  Constructs a vertex whose geometric embedding is the point `p` and pointing to face `f`.
  `Face_handle` is the type defined the triangulation data structure used in the Delaunay triangulation
  on the sphere.
  */
  TriangulationOnSphereVertexBase_2(Point p, Face_handle f);

  /// @}

  /// \name Access Functions
  /// @{

  /*!
  Returns the point.
  */
  Point point() const;

  /// @}

  /// \name Setting
  /// @{

  /*!
  Sets the point.
  */
  void set_point(Point p);

  /*!
  Inputs the non-combinatorial information given by the vertex:
  the point and other possible information.
  */
  istream& operator>>(istream& is, TriangulationOnSphereVertexBase_2& v);

  /*!
  Outputs the non combinatorial operation given by the vertex: the
  point and other possible information.
  */
  ostream& operator<<(ostream& os, const TriangulationOnSphereVertexBase_2& v);

  /// @}
};
