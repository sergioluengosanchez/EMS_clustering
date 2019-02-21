#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iterator>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Random.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <boost/lexical_cast.hpp>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_iterator vertex_iterator;
typedef Graph_traits::face_iterator face_iterator;

using namespace Rcpp;
// [[Rcpp::export]]
std::vector<double> geodesic_distance(std::string input_path, int vertex_idx)
{
  std::cout << "Reading the mesh for geodesic distance computation"<<std::endl;
  // read input polyhedron
  Triangle_mesh tmesh;
  std::ifstream input(input_path);
  input >> tmesh;
  input.close();
  std::cout << "File readed"<<std::endl;
  // initialize indices of vertices, halfedges and faces
  CGAL::set_halfedgeds_items_id(tmesh);
  const int target_vertex_index = vertex_idx;
  vertex_iterator vertex_it = vertices(tmesh).first;
  std::advance(vertex_it,vertex_idx);

  std::cout << "Computing geodesic distance "<<std::endl;
  // construct a shortest path query object and add a source point
  Surface_mesh_shortest_path shortest_paths(tmesh);
  shortest_paths.add_source_point(*vertex_it);
  // For all vertices in the tmesh, compute the points of
  // the shortest path to the source point and write them
  // into a file readable using the CGAL Polyhedron demo
  std::vector<double> distance;
  vertex_iterator vit, vit_end;
  for (boost::tie(vit, vit_end) = vertices(tmesh);
       vit != vit_end; ++vit)
  {
    distance.push_back(shortest_paths.shortest_distance_to_source_points(*vit).first);
  }
  std::cout << "Returning geodesic distance "<<std::endl;
  return distance;
}
