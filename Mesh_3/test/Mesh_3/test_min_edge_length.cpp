#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

template <typename K, typename Polyhedron>
struct Tester {
  // Domain
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K, Polyhedron> Mesh_domain;


  // Triangulation
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;

  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  int operator()(const std::string fname, const std::string out_fname)
  {
    std::ifstream input(fname);
    using namespace CGAL::parameters;

    Polyhedron polyhedron;
    input >> polyhedron;
    if (input.fail()) {
      std::cerr << "Error: Cannot read file " << fname << std::endl;
      return EXIT_FAILURE;
    }

    if (!CGAL::is_triangle_mesh(polyhedron)) {
      std::cerr << "Input geometry is not triangulated." << std::endl;
      return EXIT_FAILURE;
    }

    // Create domain
    Mesh_domain domain(polyhedron);

    domain.detect_features();

    // Mesh criteria
    Mesh_criteria criteria(edge_size = 0.2,
                           edge_min_size = 0.1,
                           facet_distance = 0.02,
                           cell_radius_edge_ratio = 3);

    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

    // Output
    CGAL::dump_c3t3(c3t3, out_fname);

    return EXIT_SUCCESS;
  }
};

int main(int argc, char* argv[])
{
  const std::string fname = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/star.off");

  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Surface_mesh<K::Point_3>                      Surface_mesh;

  return Tester<K, Surface_mesh>()(fname, "out-star");
}
