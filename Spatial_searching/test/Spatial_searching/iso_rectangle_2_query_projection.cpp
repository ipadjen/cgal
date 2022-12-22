#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_3.h>

#include <CGAL/Timer.h>
#include <CGAL/hilbert_sort.h>

// Point_3 to Point_2 projection on the fly
template<class K>
struct Projection_xy_property_map
{
  typedef typename K::Point_3 key_type;
  typedef typename K::Point_2 value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  friend value_type get(Projection_xy_property_map<K>, const key_type& k)
  {
    return value_type(k.x(), k.y());
  }
};

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;

typedef CGAL::Random_points_in_cube_3<Point_3> Random_points_iterator;
typedef CGAL::Counting_iterator<Random_points_iterator> N_Random_points_iterator;

// Search for projected points
typedef CGAL::Search_traits_2<K>Traits_base;
typedef CGAL::Search_traits_adapter<Point_3, Projection_xy_property_map<K>, Traits_base> Traits2;
typedef CGAL::Kd_tree<Traits2> Tree2;
typedef CGAL::Fuzzy_iso_box<Traits2> Fuzzy_iso_box2;

// Search directly from 3D points
typedef CGAL::Search_traits_3<K>Traits3;
typedef CGAL::Kd_tree<Traits3> Tree3;
typedef CGAL::Fuzzy_iso_box<Traits3> Fuzzy_iso_box3;

int main()
{
  CGAL::Timer timer;
  double largnum = 1e7;
  const int N = 1000000;
  const int bucket_size = 100;

  Random_points_iterator rpg;
  std::vector<Point_3> points;
    for(int i = 0; i < N; i++)
        points.push_back(*rpg++);

  //-- Test the 3D projection
  timer.start();
  Tree2 tree2(points.begin(), points.end(), bucket_size);
  timer.stop();

  // define 2D range query
  Point_2 p(0.2, 0.2);
  Point_2 q(0.7, 0.7);
  Fuzzy_iso_box2 exact_range2(p,q);

  std::vector<Point_3> result2;
  timer.start();
  tree2.search( std::back_inserter( result2 ), exact_range2);
  timer.stop();

  // sort the resulting vector for comparison
  CGAL::hilbert_sort(result2.begin(), result2.end());

//  std::cout << "The points in the box [0.2, 0.7]^2 are: " << std::endl;
//  std::copy (result.begin(), result.end(), std::ostream_iterator<Point_3>(std::cout,"\n") );
//  std::cout << std::endl;
  std::cout << "\nTime to run projected tree insert and search: " << timer.time() << " s" << std::endl;

  timer.reset();

  //-- Test directly doing 3D spatial search
  timer.start();
  Tree3 tree3(points.begin(), points.end(), bucket_size);
  timer.stop();

  // define 3D range query
  Point_3 p3(0.2, 0.2, -largnum);
  Point_3 q3(0.7, 0.7, largnum);
  Fuzzy_iso_box3 exact_range3(p3,q3);

  std::vector<Point_3> result3;
  timer.start();
  tree3.search( std::back_inserter( result3 ), exact_range3);
  timer.stop();

  // sort the resulting vector for comparison
  CGAL::hilbert_sort(result3.begin(), result3.end());

//  std::cout << "The points in the box [0.2, 0.7]^2 are: " << std::endl;
//  std::copy (result3.begin(), result3.end(), std::ostream_iterator<Point_3>(std::cout,"\n") );
//  std::cout << std::endl;
  std::cout << "\nTime to run direct tree insert and search: " << timer.time() << " s" << std::endl;

  // Compare the results between the two
  if (result2 == result3) {
      std::cout << "\nThe two resulting point vectors are equal!" << std::endl;
  } else {
      std::cout << "\nThe two resulting point vectors are NOT equal!" << std::endl;
  }

  return 0;
}
