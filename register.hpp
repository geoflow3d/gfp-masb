#include "masb_nodes.hpp"

using namespace geoflow::nodes::mat;

void register_nodes(geoflow::NodeRegister& node_register) {
    node_register.register_node<ComputeMedialAxisNode>("ComputeMedialAxisNode");
    node_register.register_node<ComputeNormalsNode>("ComputeNormalsNode");
    node_register.register_node<SegmentMakerNode>("SegmentMaker");
    node_register.register_node<TestPointsNode>("TestPoints");
    node_register.register_node<RegionGrowMedialAxisNode>("RegionGrowMedialAxis");
}