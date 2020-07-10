#define CGAL_TRACE_STREAM std::cerr

#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Octree/Walker_criterion.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>

#include <cassert>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;
typedef CGAL::Octree::Octree
        <Point_set, typename Point_set::Point_map>
        Octree;

bool test_preorder_1_node() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto nodes = octree.walk<CGAL::Octree::Walker::Preorder>();

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());

  return true;
}

bool test_preorder_9_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({-1, -1, -1});
  points.insert({1, -1, -1});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto nodes = octree.walk<CGAL::Octree::Walker::Preorder>();

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert((*iter == octree.root()[i]));
  }

  return true;
}

bool test_preorder_25_nodes() {

  // Define the dataset
  Point_set points;
  points.insert({1, 1, 1});
  points.insert({1, 1, 2});
  points.insert({1, 1, 3});
  points.insert({1, 1, 4});
  auto point_map = points.point_map();

  // Create the octree
  Octree octree(points, point_map);
  octree.refine(10, 1);

  // Create the range
  auto nodes = octree.walk<CGAL::Octree::Walker::Preorder>();

  // Check each item in the range
  auto iter = nodes.begin();
  assert(*iter == octree.root());
  iter++;
  assert((*iter == octree.root()[0]));
  iter++;
  assert((*iter == octree.root()[1]));
  iter++;
  assert((*iter == octree.root()[2]));
  iter++;
  assert((*iter == octree.root()[3]));
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert((*iter == octree.root()[3][i]));
  }
  iter++;
  assert((*iter == octree.root()[4]));
  iter++;
  assert((*iter == octree.root()[5]));
  iter++;
  assert((*iter == octree.root()[6]));
  iter++;
  assert((*iter == octree.root()[7]));
  for (int i = 0; i < 8; ++i) {
    iter++;
    assert((*iter == octree.root()[7][i]));
  }

  std::cout << octree;

  return true;
}

int main(void) {

  test_preorder_1_node();
  test_preorder_9_nodes();
  test_preorder_25_nodes();

  return 0;
}
