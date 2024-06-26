#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_2.h>
#include <CGAL/create_straight_skeleton_2.h>
#include <CGAL/Straight_skeleton_2/IO/print.h>
#include <CGAL/draw_straight_skeleton_2.h>

#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K ;

typedef K::Point_2                   Point ;
typedef CGAL::Polygon_2<K>           Polygon_2 ;
typedef CGAL::Straight_skeleton_2<K> Ss ;

typedef std::shared_ptr<Ss> SsPtr ;

int main()
{
  Polygon_2 poly ;
  poly.push_back( Point(-1,-1) ) ;
  poly.push_back( Point(0,-12) ) ;
  poly.push_back( Point(1,-1) ) ;
  poly.push_back( Point(12,0) ) ;
  poly.push_back( Point(1,1) ) ;
  poly.push_back( Point(0,12) ) ;
  poly.push_back( Point(-1,1) ) ;
  poly.push_back( Point(-12,0) ) ;

  assert(poly.is_counterclockwise_oriented());

  // You can pass the polygon via an iterator pair
  SsPtr iss = CGAL::create_interior_straight_skeleton_2(poly.vertices_begin(), poly.vertices_end());

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*iss);
  CGAL::draw(*iss);

  // Or you can pass the polygon directly, as below.

  // To create an exterior straight skeleton you need to specify a maximum offset.
  double lMaxOffset = 5 ;
  SsPtr oss = CGAL::create_exterior_straight_skeleton_2(lMaxOffset, poly);

  CGAL::Straight_skeletons_2::IO::print_straight_skeleton(*oss);
  CGAL::draw(*oss);

  return EXIT_SUCCESS;
}
