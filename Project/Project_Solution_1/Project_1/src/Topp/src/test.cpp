#include "polymap.hpp"
#include <iostream>
#include <ros/ros.h>


int main(int argc, char **argv)
{
    ros::init(argc, argv, "polymap_test");
    ros::NodeHandle nh;
    PolyMap polymap(nh, Eigen::Vector3d(10.0, 10.0, 10.0), 20);
    polymap.generate_obstacles();
    
    auto closest_pair = polymap.getSignedDistance(Eigen::Vector3d(0.0, 0.0, 0.0));
    std::cout << "Closest pair: " << closest_pair.first.transpose() << " with distance " << closest_pair.second << std::endl;
    ros::Rate lr(10);
    while (ros::ok())
    {
        polymap.visualize_obstacles();
        auto closest_pair1 = polymap.getSignedDistance(Eigen::Vector3d(0.0, 0.0, 0.0));
        ros::spinOnce();
        lr.sleep();
    }
    return 0;
}