#include "masb_nodes.hpp"
#include "region_grower_testers.hpp"

#include <cmath>
#include <algorithm>

#include <fstream>
#include <iomanip>

namespace geoflow::nodes::mat {

void ComputeMedialAxisNode::process(){
  auto point_collection = input("points").get<PointCollection>();
  auto normals_vec3f = input("normals").get<vec3f>();

  // prepare data structures and transfer data
  masb::ma_data madata;
  madata.m = point_collection.size();
  
  masb::PointList coords;
  coords.reserve(madata.m);
  for(auto& p : point_collection) {
    coords.push_back(masb::Point(p.data()));
  }
  masb::VectorList normals;
  normals.reserve(madata.m);
  for(auto& n : normals_vec3f) {
    normals.push_back(masb::Vector(n.data()));
  }
  masb::PointList ma_coords_(madata.m*2);
  std::vector<int> ma_qidx_(madata.m*2);
  
  madata.coords = &coords;
  madata.normals = &normals;
  madata.ma_coords = &ma_coords_;
  madata.ma_qidx = ma_qidx_.data();

  // compute mat points
  masb::compute_masb_points(params, madata);

  // retrieve mat points
  vec1i ma_qidx;
  ma_qidx.reserve(madata.m*2);
  for(size_t i=0 ; i<madata.m*2; ++i) {
    ma_qidx.push_back(madata.ma_qidx[i]);
  }

  PointCollection ma_coords;
  ma_coords.reserve(madata.m*2);
  for(auto& c : ma_coords_) {
    ma_coords.push_back({c[0], c[1], c[2]});
  }

  // Compute medial geometry
  vec1f ma_radii(madata.m*2);
  vec1f ma_sepangle(madata.m*2);
  vec3f ma_spoke_f1(madata.m*2);
  vec3f ma_spoke_f2(madata.m*2);
  vec3f ma_bisector(madata.m*2);
  vec3f ma_spokecross(madata.m*2);
  for(size_t i=0; i<madata.m*2; ++i) {
    auto i_ = i%madata.m;
    auto& c = ma_coords_[i];
    // feature points
    auto& f1 = coords[i_];
    auto& f2 = coords[ma_qidx[i]];
    // radius
    ma_radii[i] = Vrui::Geometry::dist(f1, c);
    // spoke vectors
    auto s1 = f1-c;
    auto s2 = f2-c;
    ma_spoke_f1[i] = {s1[0], s1[1], s1[2]};
    ma_spoke_f2[i] = {s2[0], s2[1], s2[2]};
    // bisector
    s1.normalize();
    s2.normalize();
    auto b = (s1+s2).normalize();
    ma_bisector[i] = {b[0], b[1], b[2]};
    // separation angle
    ma_sepangle[i] = std::acos(s1*s2);
    // cross product of spoke vectors
    auto scross = Vrui::Geometry::cross(s1,s2).normalize();
    ma_spokecross[i] = {scross[0], scross[1], scross[2]};
  }
  vec1i ma_is_interior(madata.m*2, 0);
  std::fill_n(ma_is_interior.begin(), madata.m, 1);

  output("ma_coords").set(ma_coords);
  output("ma_qidx").set(ma_qidx);
  output("ma_radii").set(ma_radii);
  output("ma_is_interior").set(ma_is_interior);
  output("ma_sepangle").set(ma_sepangle);
  output("ma_bisector").set(ma_bisector);
  output("ma_spoke_f1").set(ma_spoke_f1);
  output("ma_spoke_f2").set(ma_spoke_f2);
  output("ma_spokecross").set(ma_spokecross);
}


void ComputeNormalsNode::process(){
  auto point_collection = input("points").get<PointCollection>();

  masb::ma_data madata;
  madata.m = point_collection.size();
  masb::PointList coords;
  coords.reserve(madata.m);
  for(auto& p : point_collection) {
    coords.push_back(masb::Point(p.data()));
  }
  masb::VectorList normals(madata.m);
  madata.coords = &coords;
  madata.normals = &normals;

  masb::compute_normals(params, madata);

  vec3f normals_vec3f;
  normals_vec3f.reserve(madata.m);
  for(auto& n : *madata.normals) {
    normals_vec3f.push_back({n[0], n[1], n[2]});
  }

  output("normals").set(normals_vec3f);
}


void SegmentMakerNode::process(){
  auto sources = input("sources").get<PointCollection>();
  auto directions = input("directions").get<vec3f>();

  if (sources.size()!=directions.size()) {
    return;
  }

  SegmentCollection segments;
  for(size_t i=0; i<sources.size(); ++i) {
    arr3f target = {sources[i][0] + directions[i][0], sources[i][1] + directions[i][1], sources[i][2] + directions[i][2]};
    segments.push_back({sources[i], target});
  }
  output("segments").set(segments);
}

void RegionGrowMedialAxisNode::process() {
  auto ma_coords = input("ma_coords").get<PointCollection>();
  auto ma_bisector = input("ma_bisector").get<vec3f>();
  auto ma_sepangle = input("ma_sepangle").get<vec1f>();
  auto ma_radii = input("ma_radii").get<vec1f>();

  regiongrower::RegionGrower<MaData,Region> R;
  R.min_segment_count = min_count;

  MaData D(ma_coords, ma_bisector, ma_sepangle, ma_radii, k);

  switch (method) {
    case 0: {
      AngleOfVectorsTester T_bisector_angle(bisector_angle);
      R.grow_regions(D, T_bisector_angle); break;
    } case 1: {
      DiffOfAnglesTester T_separation_angle(separation_angle);
      R.grow_regions(D, T_separation_angle); break;
    } case 2: {
      BallOverlapTester T_ball_overlap(ball_overlap);
      R.grow_regions(D, T_ball_overlap); break;
    } case 3: {
      CountTester T_shape_count(shape_count);
      R.grow_regions(D, T_shape_count); break;
    } default: break;
  };
  
  vec1i segment_ids;
  for(auto& region_id : R.region_ids) {
    segment_ids.push_back(int(region_id));
  }
  output("segment_ids").set(segment_ids);
}

void MATCSVWriterNode::process()
{
  auto points = input("points").get<PointCollection>();
  auto ma_coords = input("ma_coords").get<PointCollection>();
  auto radii = input("radii").get<vec1f>();
  auto sepangle = input("sepangle").get<vec1f>();

  std::ofstream f_out(filepath);
  f_out << std::fixed << std::setprecision(2);
  f_out << "x y z x_mat1 y_mat1 z_mat1 x_mat2 y_mat2 z_mat2 radius_1 radius_2 sepangle_1 sepangle_2\n";
  for (size_t i = 0; i < points.size(); ++i)
  {
    f_out
        << points[i][0] + (*manager.data_offset)[0] << " "
        << points[i][1] + (*manager.data_offset)[1] << " "
        << points[i][2] + (*manager.data_offset)[2] << " "
        << ma_coords[i][0] + (*manager.data_offset)[0] << " "
        << ma_coords[i][1] + (*manager.data_offset)[1] << " "
        << ma_coords[i][2] + (*manager.data_offset)[2] << " "
        << ma_coords[i+points.size()][0] + (*manager.data_offset)[0] << " "
        << ma_coords[i+points.size()][1] + (*manager.data_offset)[1] << " "
        << ma_coords[i+points.size()][2] + (*manager.data_offset)[2] << " "
        << radii[i] << " "
        << radii[i+points.size()] << " "
        << sepangle[i] << " "
        << sepangle[i+points.size()] << "\n";
  }
  f_out.close();
}


void PLYWriterNode::process()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef boost::tuple<Point, int> PL;
  typedef CGAL::Nth_of_tuple_property_map<0, PL> Point_map;
  // typedef CGAL::Nth_of_tuple_property_map<1, PL> Normal_map;
  typedef CGAL::Nth_of_tuple_property_map<1, PL> Label_map;
  typedef std::vector<PL> PL_vector;

  auto points = input("points").get<PointCollection>();
  auto labels = input("labels").get<vec1i>();

  PL_vector pl_points;
  pl_points.resize(points.size());
  for (size_t i = 0; i < points.size(); ++i)
  {
    pl_points[i].get<0>() = Point(points[i][0] + (*manager.data_offset)[0], points[i][1] + (*manager.data_offset)[1], points[i][2] + (*manager.data_offset)[2]);
    pl_points[i].get<1>() = labels[i];
  }

  std::ofstream f(filepath);
  if (write_binary)
    CGAL::set_binary_mode(f); // The PLY file will be written in the binary format
  else
    f << std::fixed << std::setprecision(2);

  CGAL::write_ply_points_with_properties(f, pl_points,
                                         CGAL::make_ply_point_writer(Point_map()),
                                         //  CGAL::make_ply_normal_writer (Normal_map()),
                                         std::make_pair(Label_map(), CGAL::PLY_property<int>("segment_id")));
  f.close();
}

void PLYReaderNode::process()
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef Kernel::Vector_3 Vector;
  typedef boost::tuple<Point, Vector> PL;
  typedef CGAL::Nth_of_tuple_property_map<0, PL> Point_map;
  typedef CGAL::Nth_of_tuple_property_map<1, PL> Normal_map;
  // typedef CGAL::Nth_of_tuple_property_map<1, PL> Label_map;
  typedef std::vector<PL> PN_vector;

  PN_vector pn_points;

  std::ifstream f(filepath);

  if (!f || !CGAL::read_ply_points_with_properties(f, std::back_inserter(pn_points),
                                                   CGAL::make_ply_point_reader(Point_map()),
                                                   CGAL::make_ply_normal_reader(Normal_map())))
  {
    std::cerr << "Error: cannot read file " << filepath << std::endl;
  }
  f.close();

  PointCollection points;
  vec3f normals;
  for (auto &pn : pn_points)
  {
    auto &p = boost::get<0>(pn);
    auto &n = boost::get<1>(pn);
    points.push_back({float(p.x()), float(p.y()), float(p.z())});
    normals.push_back({float(n.x()), float(n.y()), float(n.z())});
  }
  output("points").set(points);
  output("normals").set(normals);
}

}