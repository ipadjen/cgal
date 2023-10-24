#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/tetrahedral_remeshing.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef FT(Function)(const Point&);
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Remeshing
typedef CGAL::Triangulation_3<Tr::Geom_traits,
                              Tr::Triangulation_data_structure> T3_remeshing;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr>     Mesh_criteria;
typedef Mesh_criteria::Facet_criteria Facet_criteria;
typedef Mesh_criteria::Cell_criteria  Cell_criteria;

// Sizing field
struct Spherical_sizing_field
{
  typedef ::FT FT;
  typedef Point Point_3;
  typedef Mesh_domain::Index Index;
  FT operator()(const Point_3& p, const int, const Index&) const
  {
    FT sq_d_to_origin = CGAL::squared_distance(p, Point(CGAL::ORIGIN));
    return CGAL::abs(CGAL::sqrt(sq_d_to_origin) - 0.5) / 5. + 0.025;
  }
};

// Function
FT sphere_function(const Point& p)
{
  return CGAL::squared_distance(p, Point(CGAL::ORIGIN)) - 1;
}

using namespace CGAL::parameters;

int main()
{
  Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain
                        (sphere_function, K::Sphere_3(CGAL::ORIGIN, K::FT(2)));

  // Mesh criteria
  Spherical_sizing_field size;
  Mesh_criteria criteria(facet_angle(30).facet_size(0.1).facet_distance(0.025).
                         cell_radius_edge_ratio(2).cell_size(size));

  // Mesh generation
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_exude().no_perturb());

  std::cout << "Meshing done." << std::endl;

  //Remeshing : extract triangulation
  T3_remeshing t3 = CGAL::convert_to_triangulation_3(c3t3);

  //Remeshing : coarsen
  double target_edge_length = 15.;
  CGAL::tetrahedral_isotropic_remeshing(t3, size,
      number_of_iterations(2).smooth_constrained_edges(true));

  std::cout << "Remeshing done." << std::endl;

  return 0;
}
