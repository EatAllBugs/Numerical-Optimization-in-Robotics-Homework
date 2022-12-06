#ifndef VISUALIZER_HPP
#define VISUALIZER_HPP

#include "gcopter/trajectory.hpp"
#include "gcopter/quickhull.hpp"
#include "gcopter/geo_utils.hpp"

#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>

#include <ros/ros.h>
#include <std_msgs/Float64.h>
#include <geometry_msgs/Point.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <std_msgs/ColorRGBA.h>
#include <string.h>

// Visualizer for the planner
class Visualizer
{
private:
    // config contains the scale for some markers
    ros::NodeHandle nh;

    // These are publishers for path, waypoints on the trajectory,
    // the entire trajectory, the mesh of free-space polytopes,
    // the edge of free-space polytopes, and spheres for safety radius
    ros::Publisher routePub;
    ros::Publisher wayPointsPub;
    ros::Publisher trajectoryPub;
    ros::Publisher meshPub;
    ros::Publisher edgePub;
    ros::Publisher spherePub;
    ros::Publisher vistrajPub;

public:
    ros::Publisher speedPub;
    ros::Publisher thrPub;
    ros::Publisher tiltPub;
    ros::Publisher bdrPub;

    //Joeyyu: pub the exact pitch and roll angle in Z-Y-X order

    ros::Publisher pitchPub;
    ros::Publisher safePub;
    ros::Publisher rollPub;
    ros::Publisher fovPub;
    ros::Publisher dronePub;

    ros::Publisher minco_test_pub;

    //Joeyyu: Fov visuailize
    double max_dis_ = 4.0;
    // double x_max_dis_gain_ = 0.64;
    // double y_max_dis_gain_ = 0.82;
    double x_max_dis_gain_ = 0.60;
    double y_max_dis_gain_ = 0.60;
    visualization_msgs::Marker markerNode_fov;
    visualization_msgs::Marker markerEdge_fov;
    visualization_msgs::Marker marker_line, fast_marker_line;
    std::vector<Eigen::Vector3d> fov_node;

    //Joeyyu: FSC250 visualization
    visualization_msgs::Marker meshROS;
    //std::string mesh_resource = std::string("package://gcopter/meshes/swarm_drone.dae");
    double mesh_scale, relCostTol;



public:
    Visualizer(ros::NodeHandle &nh_)
        : nh(nh_)
    {
        

        routePub = nh.advertise<visualization_msgs::Marker>("/visualizer/route", 10);
        wayPointsPub = nh.advertise<visualization_msgs::Marker>("/visualizer/waypoints", 10);
        trajectoryPub = nh.advertise<visualization_msgs::Marker>("/visualizer/trajectory", 10);
        vistrajPub = nh.advertise<visualization_msgs::Marker>("/visualizer/vistraj", 10);
        meshPub = nh.advertise<visualization_msgs::Marker>("/visualizer/mesh", 1000);
        edgePub = nh.advertise<visualization_msgs::Marker>("/visualizer/edge", 1000);
        spherePub = nh.advertise<visualization_msgs::Marker>("/visualizer/spheres", 1000);
        speedPub = nh.advertise<std_msgs::Float64>("/visualizer/speed", 1000);
        thrPub = nh.advertise<std_msgs::Float64>("/visualizer/total_thrust", 1000);
        tiltPub = nh.advertise<std_msgs::Float64>("/visualizer/tilt_angle", 1000);
        bdrPub = nh.advertise<std_msgs::Float64>("/visualizer/body_rate", 1000);

        dronePub = nh.advertise<visualization_msgs::Marker>("visualizer/mesh_drone",100);
        fovPub = nh.advertise<visualization_msgs::MarkerArray>("visualizer/fov_visual",100);

        minco_test_pub = nh.advertise<visualization_msgs::MarkerArray>("visualizer/minco_test",100);


    //Joeyyu:pub the exact pitch and roll angle in Z-Y-X order
        safePub = nh.advertise<std_msgs::Float64>("/visualizer/safe_flag", 1000);
        pitchPub = nh.advertise<std_msgs::Float64>("/visualizer/pitch_angle",1000);
        rollPub = nh.advertise<std_msgs::Float64>("/visualizer/roll_angle",1000);
 
    //Joeyyu: fov initialization
        fov_visual_init("odom");
    }

    template <int D>
    inline void minco_test_trajs_pub(std::vector<Trajectory<D> > trajs)
    {
        visualization_msgs::MarkerArray trajArray;
        for (int i = 0; i< trajs.size(); i++)
        {
            visualization_msgs::Marker trajMarker;
            trajMarker.header.frame_id = "odom";
            trajMarker.header.stamp = ros::Time::now();
            trajMarker.ns = "traj";
            trajMarker.id = i;
            trajMarker.type = visualization_msgs::Marker::LINE_STRIP;
            trajMarker.action = visualization_msgs::Marker::ADD;
            trajMarker.pose.orientation.w = 1.0;
            trajMarker.scale.x = 0.1;
            trajMarker.scale.y = 0.1;
            trajMarker.scale.z = 0.1;
            trajMarker.color.a = 1.0;
            trajMarker.color.r = 0.0;
            trajMarker.color.g = 1.0;
            trajMarker.color.b = 0.0;
            trajMarker.lifetime = ros::Duration(0);
            trajMarker.frame_locked = false;
            trajMarker.points.clear();

            Trajectory<D> traj = trajs[i];
            if (traj.getPieceNum() > 0)
            {
                double vel_max = 10.0;
                double T = 0.01;
                Eigen::Vector3d lastX = traj.getPos(0.0);
                Eigen::Vector3d lastVel = traj.getVel(0.0);
                std_msgs::ColorRGBA color_c, color_last;
                for (double t = T; t < traj.getTotalDuration(); t += T)
                {
                    geometry_msgs::Point point;
                    Eigen::Vector3d X = traj.getPos(t);
                    Eigen::Vector3d Vel = traj.getVel(t);
                    point.x = lastX(0);
                    point.y = lastX(1);
                    point.z = lastX(2);
                    color_last = getColor(lastVel.norm()/vel_max);
                    trajMarker.colors.push_back(color_last);
                    trajMarker.points.push_back(point);
                    point.x = X(0);
                    point.y = X(1);
                    point.z = X(2);
                    color_c = getColor(Vel.norm()/vel_max);
                    trajMarker.colors.push_back(color_c);
                    trajMarker.points.push_back(point);
                    lastX = X;
                    lastVel = Vel;
                }
            }
            trajArray.markers.push_back(trajMarker);
        }
        minco_test_pub.publish(trajArray);


    }

    //Joeyyu: Visualize the visiable trajectory
    template <int D>
    inline void vistraj_pub(const Trajectory<D> &traj, const double &progress_t, const double & delta, const bool & safe)
    {
        visualization_msgs::Marker routeMarker, trajMarker, vistrajDeleter;

        routeMarker.id = 0;
        routeMarker.type = visualization_msgs::Marker::LINE_LIST;
        routeMarker.header.stamp = ros::Time::now();
        routeMarker.header.frame_id = "odom";
        routeMarker.pose.orientation.w = 1.00;
        routeMarker.action = visualization_msgs::Marker::ADD;
        routeMarker.ns = "route";
        routeMarker.color.r = 1.00;
        routeMarker.color.g = 0.00;
        routeMarker.color.b = 0.00;
        routeMarker.color.a = 1.00;
        routeMarker.scale.x = 0.1;

        trajMarker = routeMarker;
        trajMarker.header.frame_id = "odom";
        trajMarker.id = 0;
        trajMarker.ns = "vistraj";
        trajMarker.color.r = safe? 1.00: 0.00;
        trajMarker.color.g = 0.00;
        trajMarker.color.b = 0.00;
        trajMarker.scale.x = 0.50;

        vistrajDeleter = trajMarker;
        vistrajDeleter.action = visualization_msgs::Marker::DELETEALL;

        if (traj.getPieceNum() > 0)
        {
            double T = 0.01;
            Eigen::Vector3d lastX = traj.getPos(delta);
            for (double t = delta+T; t <= progress_t; t += T)
            {
                geometry_msgs::Point point;
                Eigen::Vector3d X = traj.getPos(t);
                point.x = lastX(0);
                point.y = lastX(1);
                point.z = lastX(2);
                trajMarker.points.push_back(point);
                point.x = X(0);
                point.y = X(1);
                point.z = X(2);
                trajMarker.points.push_back(point);
                lastX = X;
            }
            vistrajPub.publish(vistrajDeleter);
            vistrajPub.publish(trajMarker);
        }






    } 

    // Visualize the trajectory and its front-end path
    template <int D>
    inline void visualize(const Trajectory<D> &traj,
                          const std::vector<Eigen::Vector3d> &route)
    {
        visualization_msgs::Marker routeMarker, wayPointsMarker, trajMarker;

        routeMarker.id = 0;
        routeMarker.type = visualization_msgs::Marker::LINE_LIST;
        routeMarker.header.stamp = ros::Time::now();
        routeMarker.header.frame_id = "odom";
        routeMarker.pose.orientation.w = 1.00;
        routeMarker.action = visualization_msgs::Marker::ADD;
        routeMarker.ns = "route";
        routeMarker.color.r = 1.00;
        routeMarker.color.g = 0.00;
        routeMarker.color.b = 0.00;
        routeMarker.color.a = 1.00;
        routeMarker.scale.x = 0.1;

        wayPointsMarker = routeMarker;
        wayPointsMarker.id = -wayPointsMarker.id - 1;
        wayPointsMarker.type = visualization_msgs::Marker::SPHERE_LIST;
        wayPointsMarker.ns = "waypoints";
        wayPointsMarker.color.r = 1.00;
        wayPointsMarker.color.g = 0.00;
        wayPointsMarker.color.b = 0.00;
        wayPointsMarker.scale.x = 0.35;
        wayPointsMarker.scale.y = 0.35;
        wayPointsMarker.scale.z = 0.35;

        trajMarker = routeMarker;
        trajMarker.header.frame_id = "odom";
        trajMarker.id = 0;
        trajMarker.ns = "trajectory";
        trajMarker.color.r = 0.00;
        trajMarker.color.g = 0.50;
        trajMarker.color.b = 1.00;
        trajMarker.scale.x = 0.30;

        if (route.size() > 0)
        {
            bool first = true;
            Eigen::Vector3d last;
            for (auto it : route)
            {
                if (first)
                {
                    first = false;
                    last = it;
                    continue;
                }
                geometry_msgs::Point point;

                point.x = last(0);
                point.y = last(1);
                point.z = last(2);
                routeMarker.points.push_back(point);
                point.x = it(0);
                point.y = it(1);
                point.z = it(2);
                routeMarker.points.push_back(point);
                last = it;
            }

            routePub.publish(routeMarker);
        }

        if (traj.getPieceNum() > 0)
        {
            Eigen::MatrixXd wps = traj.getPositions();
            for (int i = 0; i < wps.cols(); i++)
            {
                geometry_msgs::Point point;
                point.x = wps.col(i)(0);
                point.y = wps.col(i)(1);
                point.z = wps.col(i)(2);
                wayPointsMarker.points.push_back(point);
            }

            wayPointsPub.publish(wayPointsMarker);
        }

        if (traj.getPieceNum() > 0)
        {
            double vel_max = 10.0;
            double T = 0.01;
            Eigen::Vector3d lastX = traj.getPos(0.0);
            Eigen::Vector3d lastVel = traj.getVel(0.0);
            std_msgs::ColorRGBA color_c, color_last;
            for (double t = T; t < traj.getTotalDuration(); t += T)
            {
                geometry_msgs::Point point;
                Eigen::Vector3d X = traj.getPos(t);
                Eigen::Vector3d Vel = traj.getVel(t);
                point.x = lastX(0);
                point.y = lastX(1);
                point.z = lastX(2);
                color_last = getColor(lastVel.norm()/vel_max);
                trajMarker.colors.push_back(color_last);
                trajMarker.points.push_back(point);
                point.x = X(0);
                point.y = X(1);
                point.z = X(2);
                color_c = getColor(Vel.norm()/vel_max);
                trajMarker.colors.push_back(color_c);
                trajMarker.points.push_back(point);
                lastX = X;
                lastVel = Vel;
            }
            trajectoryPub.publish(trajMarker);
        }
    }

    //Joeyyu: fov visualization
    inline void fov_visual_init(std::string msg_frame_id) 
    {
        double x_max_dis = max_dis_ * x_max_dis_gain_;
        double y_max_dis = max_dis_ * y_max_dis_gain_;

        fov_node.resize(5);
        fov_node[0][0] = 0;
        fov_node[0][1] = 0;
        fov_node[0][2] = 0;

        fov_node[1][2] = x_max_dis;
        fov_node[1][1] = y_max_dis;
        fov_node[1][0] = max_dis_;

        fov_node[2][2] = x_max_dis;
        fov_node[2][1] = -y_max_dis;
        fov_node[2][0] = max_dis_;

        fov_node[3][2] = -x_max_dis;
        fov_node[3][1] = -y_max_dis;
        fov_node[3][0] = max_dis_;

        fov_node[4][2] = -x_max_dis;
        fov_node[4][1] = y_max_dis;
        fov_node[4][0] = max_dis_;

        markerNode_fov.header.frame_id = msg_frame_id;
        // markerNode_fov.header.stamp = msg_time;
        markerNode_fov.action = visualization_msgs::Marker::ADD;
        markerNode_fov.type = visualization_msgs::Marker::SPHERE_LIST;
        markerNode_fov.ns = "fov_nodes";
        // markerNode_fov.id = 0;
        markerNode_fov.pose.orientation.w = 1;
        markerNode_fov.scale.x = 0.05;
        markerNode_fov.scale.y = 0.05;
        markerNode_fov.scale.z = 0.05;
        markerNode_fov.color.r = 0;
        markerNode_fov.color.g = 0.8;
        markerNode_fov.color.b = 1;
        markerNode_fov.color.a = 1;

        markerEdge_fov.header.frame_id = msg_frame_id;
        // markerEdge_fov.header.stamp = msg_time;
        markerEdge_fov.action = visualization_msgs::Marker::ADD;
        markerEdge_fov.type = visualization_msgs::Marker::LINE_LIST;
        markerEdge_fov.ns = "fov_edges";
        // markerEdge_fov.id = 0;
        markerEdge_fov.pose.orientation.w = 1;
        markerEdge_fov.scale.x = 0.05;
        markerEdge_fov.color.r = 0.5f;
        markerEdge_fov.color.g = 0.0;
        markerEdge_fov.color.b = 0.0;
        markerEdge_fov.color.a = 1;
    }

    inline void pub_fov_visual(const Eigen::Vector3d& p, Eigen::Quaterniond& q) 
    {

        visualization_msgs::Marker clear_previous_msg;
        clear_previous_msg.action = visualization_msgs::Marker::DELETEALL;

        visualization_msgs::MarkerArray markerArray_fov;
        markerNode_fov.points.clear();
        markerEdge_fov.points.clear();

        std::vector<geometry_msgs::Point> fov_node_marker;
        for (int i = 0; i < (int)fov_node.size(); i++) {
            Eigen::Vector3d vector_temp;
            vector_temp = q * fov_node[i] + p;
            geometry_msgs::Point point_temp;
            point_temp.x = vector_temp[0];
            point_temp.y = vector_temp[1];
            point_temp.z = vector_temp[2];
            fov_node_marker.push_back(point_temp);
        }


        markerNode_fov.points.push_back(fov_node_marker[0]);
        markerNode_fov.points.push_back(fov_node_marker[1]);
        markerNode_fov.points.push_back(fov_node_marker[2]);
        markerNode_fov.points.push_back(fov_node_marker[3]);
        markerNode_fov.points.push_back(fov_node_marker[4]);

        markerEdge_fov.points.push_back(fov_node_marker[0]);
        markerEdge_fov.points.push_back(fov_node_marker[1]);

        markerEdge_fov.points.push_back(fov_node_marker[0]);
        markerEdge_fov.points.push_back(fov_node_marker[2]);

        markerEdge_fov.points.push_back(fov_node_marker[0]);
        markerEdge_fov.points.push_back(fov_node_marker[3]);

        markerEdge_fov.points.push_back(fov_node_marker[0]);
        markerEdge_fov.points.push_back(fov_node_marker[4]);

        markerEdge_fov.points.push_back(fov_node_marker[0]);
        markerEdge_fov.points.push_back(fov_node_marker[1]);

        markerEdge_fov.points.push_back(fov_node_marker[1]);
        markerEdge_fov.points.push_back(fov_node_marker[2]);

        markerEdge_fov.points.push_back(fov_node_marker[2]);
        markerEdge_fov.points.push_back(fov_node_marker[3]);

        markerEdge_fov.points.push_back(fov_node_marker[3]);
        markerEdge_fov.points.push_back(fov_node_marker[4]);

        markerEdge_fov.points.push_back(fov_node_marker[4]);
        markerEdge_fov.points.push_back(fov_node_marker[1]);

        markerArray_fov.markers.push_back(clear_previous_msg);
        markerArray_fov.markers.push_back(markerNode_fov);
        markerArray_fov.markers.push_back(markerEdge_fov);
        fovPub.publish(markerArray_fov);
    }

    inline void pub_mesh_drone(const Eigen::Vector3d& pose, Eigen::Vector4d& q, double scale, std::string mesh_resource)
    {
        
        meshROS.header.frame_id = "odom";
        meshROS.header.stamp = ros::Time::now();
        meshROS.ns = "mesh";
        meshROS.id = 0;
        meshROS.type = visualization_msgs::Marker::MESH_RESOURCE;
        meshROS.action = visualization_msgs::Marker::ADD;
        meshROS.mesh_use_embedded_materials = true;

        meshROS.pose.position.x = pose(0);
        meshROS.pose.position.y = pose(1);
        meshROS.pose.position.z = pose(2);

        meshROS.pose.orientation.w = q(0);
        meshROS.pose.orientation.x = q(1);
        meshROS.pose.orientation.y = q(2);
        meshROS.pose.orientation.z = q(3);
        meshROS.scale.x = scale;
        meshROS.scale.y = scale;
        meshROS.scale.z = scale;

        meshROS.color.a = 0;
        meshROS.color.r = 0;
        meshROS.color.g = 0;
        meshROS.color.b = 0;
        meshROS.mesh_resource = mesh_resource;
        dronePub.publish(meshROS);




    }



    // Visualize some polytopes in H-representation
    inline void visualizePolytope(const std::vector<Eigen::MatrixX4d> &hPolys)
    {

        // Due to the fact that H-representation cannot be directly visualized
        // We first conduct vertex enumeration of them, then apply quickhull
        // to obtain triangle meshs of polyhedra
        Eigen::Matrix3Xd mesh(3, 0), curTris(3, 0), oldTris(3, 0);
        for (size_t id = 0; id < hPolys.size(); id++)
        {
            oldTris = mesh;
            Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vPoly;
            geo_utils::enumerateVs(hPolys[id], vPoly);

            quickhull::QuickHull<double> tinyQH;
            const auto polyHull = tinyQH.getConvexHull(vPoly.data(), vPoly.cols(), false, true);
            const auto &idxBuffer = polyHull.getIndexBuffer();
            int hNum = idxBuffer.size() / 3;

            curTris.resize(3, hNum * 3);
            for (int i = 0; i < hNum * 3; i++)
            {
                curTris.col(i) = vPoly.col(idxBuffer[i]);
            }
            mesh.resize(3, oldTris.cols() + curTris.cols());
            mesh.leftCols(oldTris.cols()) = oldTris;
            mesh.rightCols(curTris.cols()) = curTris;
        }

        // RVIZ support tris for visualization
        visualization_msgs::Marker meshMarker, edgeMarker;

        meshMarker.id = 0;
        meshMarker.header.stamp = ros::Time::now();
        meshMarker.header.frame_id = "odom";
        meshMarker.pose.orientation.w = 1.00;
        meshMarker.action = visualization_msgs::Marker::ADD;
        meshMarker.type = visualization_msgs::Marker::TRIANGLE_LIST;
        meshMarker.ns = "mesh";
        meshMarker.color.r = 0.00;
        meshMarker.color.g = 0.00;
        meshMarker.color.b = 1.00;
        meshMarker.color.a = 0.15;
        meshMarker.scale.x = 1.0;
        meshMarker.scale.y = 1.0;
        meshMarker.scale.z = 1.0;

        edgeMarker = meshMarker;
        edgeMarker.type = visualization_msgs::Marker::LINE_LIST;
        edgeMarker.ns = "edge";
        edgeMarker.color.r = 0.00;
        edgeMarker.color.g = 1.00;
        edgeMarker.color.b = 1.00;
        edgeMarker.color.a = 1.00;
        edgeMarker.scale.x = 0.02;

        geometry_msgs::Point point;

        int ptnum = mesh.cols();

        for (int i = 0; i < ptnum; i++)
        {
            point.x = mesh(0, i);
            point.y = mesh(1, i);
            point.z = mesh(2, i);
            meshMarker.points.push_back(point);
        }

        for (int i = 0; i < ptnum / 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                point.x = mesh(0, 3 * i + j);
                point.y = mesh(1, 3 * i + j);
                point.z = mesh(2, 3 * i + j);
                edgeMarker.points.push_back(point);
                point.x = mesh(0, 3 * i + (j + 1) % 3);
                point.y = mesh(1, 3 * i + (j + 1) % 3);
                point.z = mesh(2, 3 * i + (j + 1) % 3);
                edgeMarker.points.push_back(point);
            }
        }

        meshPub.publish(meshMarker);
        edgePub.publish(edgeMarker);

        return;
    }




    // Visualize all spheres with centers sphs and the same radius
    inline void visualizeSphere(const Eigen::Vector3d &center,
                                const double &radius)
    {
        visualization_msgs::Marker sphereMarkers, sphereDeleter;

        sphereMarkers.id = 0;
        sphereMarkers.type = visualization_msgs::Marker::SPHERE_LIST;
        sphereMarkers.header.stamp = ros::Time::now();
        sphereMarkers.header.frame_id = "odom";
        sphereMarkers.pose.orientation.w = 1.00;
        sphereMarkers.action = visualization_msgs::Marker::ADD;
        sphereMarkers.ns = "spheres";
        sphereMarkers.color.r = 0.00;
        sphereMarkers.color.g = 0.00;
        sphereMarkers.color.b = 1.00;
        sphereMarkers.color.a = 1.00;
        sphereMarkers.scale.x = radius * 2.0;
        sphereMarkers.scale.y = radius * 2.0;
        sphereMarkers.scale.z = radius * 2.0;

        sphereDeleter = sphereMarkers;
        sphereDeleter.action = visualization_msgs::Marker::DELETE;

        geometry_msgs::Point point;
        point.x = center(0);
        point.y = center(1);
        point.z = center(2);
        sphereMarkers.points.push_back(point);

        //spherePub.publish(sphereDeleter);
        //spherePub.publish(sphereMarkers);
    }

    inline void visualizeStartGoal(const Eigen::Vector3d &center,
                                   const double &radius,
                                   const int sg)
    {
        visualization_msgs::Marker sphereMarkers, sphereDeleter;

        sphereMarkers.id = sg;
        sphereMarkers.type = visualization_msgs::Marker::SPHERE_LIST;
        sphereMarkers.header.stamp = ros::Time::now();
        sphereMarkers.header.frame_id = "odom";
        sphereMarkers.pose.orientation.w = 1.00;
        sphereMarkers.action = visualization_msgs::Marker::ADD;
        sphereMarkers.ns = "StartGoal";
        sphereMarkers.color.r = 1.00;
        sphereMarkers.color.g = 0.00;
        sphereMarkers.color.b = 0.00;
        sphereMarkers.color.a = 1.00;
        sphereMarkers.scale.x = radius * 2.0;
        sphereMarkers.scale.y = radius * 2.0;
        sphereMarkers.scale.z = radius * 2.0;

        sphereDeleter = sphereMarkers;
        sphereDeleter.action = visualization_msgs::Marker::DELETEALL;

        geometry_msgs::Point point;
        point.x = center(0);
        point.y = center(1);
        point.z = center(2);
        sphereMarkers.points.push_back(point);

        if (sg == 0)
        {
            spherePub.publish(sphereDeleter);
            ros::Duration(1.0e-9).sleep();
            sphereMarkers.header.stamp = ros::Time::now();
        }
        spherePub.publish(sphereMarkers);
    }

    inline std_msgs::ColorRGBA getColor(double portial)
    {
        if(portial > 1) portial = 1;
        double r,g,b,a;

        if(portial < 0.25)
        {
            r = 0.0;
            g = 0.0 + portial/0.25;
            b = 1.0; 
            
        }
        else if(0.25 <= portial < 0.5)
        {   
            r = 0.0;
            g = 1.0;
            b = 1.0 - (portial-0.25)/0.25;

        }
        else if(0.5 <= portial < 0.75)
        {
            r = 0.0 + (portial-0.5)/0.25;
            g = 1.0;
            b = 0.0;
        }
        else
        {
            r = 1.0;
            g = 1.0 - (portial-0.75)/0.25;
            b = 0.0;
        }

        std_msgs::ColorRGBA ret;
        ret.r = r;
        ret.g = g;
        ret.b = b;
        return ret;
    }
};

#endif