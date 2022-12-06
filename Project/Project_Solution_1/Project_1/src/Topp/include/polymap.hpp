#ifndef POLYMAP_HPP_
#define POLYMAP_HPP_

#include <Eigen/Eigen>
#include <vector>
#include <iostream>
#include <random>
#include <cstdlib>
#include <climits>
#include "support/quickhull.hpp"
#include "support/geo_utils.hpp"
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <geometry_msgs/Point.h>
#include "support/sdqp.hpp""



class PolyMap 
{
private:
    /* data */
    int num_obstacles;
    Eigen::Vector3d map_size;
    std::vector<Eigen::MatrixX4d> hPolys;

    ros::NodeHandle nh;
    ros::Publisher meshPub;
    ros::Publisher edgePub;
    ros::Publisher disvectorPub;

    

public:
    PolyMap(ros::NodeHandle &nh_, Eigen::Vector3d mp_size, int number_of_obstacles):nh(nh_)
    {
        map_size = mp_size;
        num_obstacles = number_of_obstacles;
        //hPolys.resize(num_obstacles);
        meshPub = nh.advertise<visualization_msgs::Marker>("polymap_mesh", 1);
        edgePub = nh.advertise<visualization_msgs::Marker>("polymap_edges", 1);
        disvectorPub = nh.advertise<visualization_msgs::Marker>("polymap_disvector", 1);
        

    }

    inline void generate_obstacles();
    inline std::vector<Eigen::MatrixX4d> get_obstacles(){return hPolys;}
    inline void visualize_obstacles();
    inline std::pair<Eigen::Vector3d, double> getSignedDistance(Eigen::Vector3d point);
    inline void getPotentialCostGrad(const Eigen::VectorXd xi, Eigen::VectorXd &g, double &f);

    
};

    inline void PolyMap::generate_obstacles()
    {
        for(int i = 0; i<num_obstacles; i++)
        {
          
        

        Eigen::Matrix3Xd mesh, recoveredV;
        Eigen::Matrix<double, 3, -1, Eigen::ColMajor> vertices;
        Eigen::Vector3d inner;

        // Randomly generate a set of points
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-0.5, 1.0);
        std::uniform_real_distribution<> dis2(-1.0, 1.0);
        vertices.resize(3, 5);
        for (int i = 0; i < 5; i++)
        {
            vertices.col(i) << dis2(gen), dis2(gen), dis(gen);
            //std::cout << vertices.col(i).transpose() << std::endl;
        }
        vertices.array() *=1.5;
        vertices.row(0).array() += dis2(gen) * 5.0 * 1.0;
        vertices.row(1).array() += dis2(gen) * 5.0 * 1.0;
        vertices.row(2).array() += dis(gen) * 2.0 * 1.0;

        // Get the convex hull in mesh
        quickhull::QuickHull<double> qh;
        const auto cvxHull = qh.getConvexHull(vertices.data(),
                                            vertices.cols(),
                                            false, false);
        const auto &idBuffer = cvxHull.getIndexBuffer();
        const auto &vtBuffer = cvxHull.getVertexBuffer();
        int ids = idBuffer.size();
        mesh.resize(3, ids);
        quickhull::Vector3<double> v;
        for (int i = 0; i < ids; i++)
        {
            v = vtBuffer[idBuffer[i]];
            mesh(0, i) = v.x;
            mesh(1, i) = v.y;
            mesh(2, i) = v.z;
        }

        // Obtain the half space intersection form from the mesh
        Eigen::MatrixX4d hPoly(ids / 3, 4);
        Eigen::Vector3d normal, point, edge0, edge1;
        for (int i = 0; i < ids / 3; i++)
        {
            point = mesh.col(3 * i + 1);
            edge0 = point - mesh.col(3 * i);
            edge1 = mesh.col(3 * i + 2) - point;
            normal = edge0.cross(edge1).normalized();
            hPoly(i, 0) = normal(0);
            hPoly(i, 1) = normal(1);
            hPoly(i, 2) = normal(2);
            hPoly(i, 3) = -normal.dot(point);
        }
        hPolys.push_back(hPoly);
        }   
        
        // Eigen::MatrixX4d hPoly2(6, 4);
        // hPoly2.row(0) << 1, 0, 0, -1;
        // hPoly2.row(1) << -1, 0, 0, -1;
        // hPoly2.row(2) << 0, 1, 0, -1;
        // hPoly2.row(3) << 0, -1, 0, -1;
        // hPoly2.row(4) << 0, 0, 1, -1;
        // hPoly2.row(5) << 0, 0, -1, -1;

        // hPolys.push_back(hPoly2);
        //}

    
    }

    inline void PolyMap::visualize_obstacles()
    {
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
        meshMarker.header.frame_id = "map";
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

    inline std::pair<Eigen::Vector3d, double> PolyMap::getSignedDistance(Eigen::Vector3d point)
    {
        Eigen::Vector3d x_final, x_tmp;
        double dist_final = 1e+10, dist_tmp_square;
        for(int i = 0; i<hPolys.size(); i++)
        {
            int hNum = hPolys[i].rows();
            Eigen::MatrixX3d A  = hPolys[i].block(0, 0, hNum, 3);
            Eigen::VectorXd b = -hPolys[i].rightCols(1);
            Eigen::Matrix3d Q = 2*Eigen::Matrix3d::Identity();
            Eigen::Vector3d c = -2*point;
            dist_tmp_square  = sdqp::sdqp(Q,c,A,b,x_tmp);
            //std::cout << "dist_tmp_square: " << dist_tmp_square << std::endl;

            if (dist_tmp_square == 0)
            {
                //std::cout<<"polymap: inside obstacle"<<std::endl;
                double dist_tmp_inside = 1e+10;
                Eigen::Vector3d x_tmp_inside;
                for (int j = 0; j < hNum; j++)
                {
                    Eigen::Vector3d aj = A.row(j).transpose();
                    double bj = b(j);
                    Eigen::Vector3d p_v = (aj.dot(point)-bj)/(aj.dot(aj)) * aj;
                    //std::cout <<"p_v.norm() = "<<p_v.norm()<<std::endl;
                    if (p_v.norm() < dist_tmp_inside)
                    {
                        dist_tmp_inside = p_v.norm();
                        x_tmp_inside = point-p_v;
                    }
                    
                }
                x_final = x_tmp_inside;
                dist_final = -dist_tmp_inside;
                break;

                //return std::make_pair(x_tmp_inside, -dist_tmp_inside);
            }
            if((x_tmp-point).norm() < dist_final)
            {
                dist_final = (x_tmp-point).norm() ;
                x_final = x_tmp;
            }

            
        }
        // visualization_msgs::Marker disvectorMarker;
        // disvectorMarker.id = 0;
        // disvectorMarker.header.stamp = ros::Time::now();
        // disvectorMarker.header.frame_id = "map";
        // disvectorMarker.action = visualization_msgs::Marker::ADD;
        // disvectorMarker.type = visualization_msgs::Marker::ARROW;
        // disvectorMarker.ns = "mesh";
        // disvectorMarker.color.r = 1.00;
        // disvectorMarker.color.g = 0.00;
        // disvectorMarker.color.b = 0.00;
        // disvectorMarker.color.a = 1.00;
        // disvectorMarker.scale.x = 0.15;
        // disvectorMarker.scale.y = 0.2;
        // geometry_msgs::Point point_start, point_end;
        // point_start.x = point(0);
        // point_start.y = point(1);
        // point_start.z = point(2);
        // point_end.x = x_final(0);
        // point_end.y = x_final(1);
        // point_end.z = x_final(2);
        // disvectorMarker.points.push_back(point_start);
        // disvectorMarker.points.push_back(point_end);
        // disvectorPub.publish(disvectorMarker);

        return std::make_pair(x_final, dist_final);
    }

    inline void PolyMap::getPotentialCostGrad(const Eigen::VectorXd xi, Eigen::VectorXd &g, double &f)
    {
    int size_n = xi.size() /3;
    Eigen::VectorXd g_tmp  = Eigen::VectorXd::Zero(xi.size());

    //int size_o = obstacle_pair.size();


    for (int i = 0; i<size_n;i++)
    {
        Eigen::Vector3d point_i(xi(3*i), xi(3*i+1), xi(3*i+2));
        //std::cout << "point_i in out_loop " << point_i.transpose() <<std::endl;
        Eigen::MatrixXd grad_xi_wrt_x = Eigen::MatrixXd::Zero(3,xi.size());
        grad_xi_wrt_x.block(0, 3*i, 3, 3) = Eigen::Matrix3d::Identity();

        auto dist_pair = getSignedDistance(point_i);
        Eigen::Vector3d point_i_dist = dist_pair.first;
        double dist_i = dist_pair.second;
        double criterion = 0.3 - dist_i;
        if (criterion >= 0)
        {
            f+= 1000* criterion;
            Eigen::Vector3d gi_wrt_xi;
            Eigen::VectorXd gi_wrt_x(xi.size());
            if (dist_i <= 0)
            {
                gi_wrt_xi = 1000 * 1/((point_i - point_i_dist).norm()) * (point_i - point_i_dist);

            }
            else
            {
                gi_wrt_xi = -1000 * 1/((point_i - point_i_dist).norm()) * (point_i - point_i_dist);
            }
            gi_wrt_x = grad_xi_wrt_x.transpose() * gi_wrt_xi;
            g_tmp   += gi_wrt_x;
           
        }

        // for (int j = 0; j < size_o; j++)
        // {
        //     //std::cout << "size_o: "  << size_o<< std::endl;
        //     Eigen::Vector3d obs_j = obstacle_pair[j].first;
        //     double r_j = obstacle_pair[j].second;

        //     //std::cout <<"obj111: " << obs_j.transpose() << std::endl;
        //     //std::cout << "r_j111: " <<  r_j << std::endl;


        //     double criterion_ij = r_j - (point_i - obs_j).norm();
        //     double smooth_factor = 0.01;
        //     double f_smooth = 0.0, g_smooth = 0.0;

        //     //if (smoothedL1(criterion_ij, smooth_factor, f_smooth, g_smooth ))
        //     if(criterion_ij >= 0)
        //     {
        //         if (print_flag)
        //         {
        //             if (criterion_ij ==0)
        //             {
        //                 std::cout <<"jjjjjjjjjjjjjjjjjjjjjj" << std::endl;
        //             } 
        //         std::cout <<"obj: " << obs_j.transpose() << std::endl;
        //         std::cout << "r_j: " <<  r_j << std::endl;
        //         std::cout << "point_i: " << point_i.transpose() << std::endl;
        //         std::cout << "criterion_ij:" << criterion_ij << std::endl;

        //         }

        //         //f += 1000 * f_smooth;
        //         f += 1000 * criterion_ij;
        //         Eigen::VectorXd g_ij_wrt_x(xi.size());
        //         Eigen::Vector3d g_ij_wrt_xi;
        //         if((point_i - obs_j).norm() == 0.0)
        //         {
        //             std::cout << "norm error================jjjjjjjjjjjjjjjjjjj===============" << std::endl;

        //         }
        //         g_ij_wrt_xi = -1000 *  0.5 * (1/ (point_i -obs_j).norm()) * (2 * point_i - 2* obs_j);
        //         g_ij_wrt_x = grad_xi_wrt_x.transpose() * g_ij_wrt_xi;

        //         //g_tmp +=  g_smooth * g_ij_wrt_x;
        //         g_tmp +=   g_ij_wrt_x;
                
        //     }
        // }
        

    }
    g += g_tmp;
}



































#endif /* POLYMAP_HPP_ */