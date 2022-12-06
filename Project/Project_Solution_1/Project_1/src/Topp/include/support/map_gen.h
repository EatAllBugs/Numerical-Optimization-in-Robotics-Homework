#ifndef MAP_GEN_H_
#define MAP_GEN_H_

#include <Eigen/Eigen>
#include <iostream>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <ros/ros.h>
#include <utility>
#include <ctime>
#include <cstdlib>




class map_gen
{
private:
    Eigen::Vector3d map_size;
    int num_obstacle;
    ros::NodeHandle nh;
    ros::Publisher obstacle_pub;
    std::vector<std::pair<Eigen::Vector3d, double>> obstacle_pair;
    bool print_flag = false;

public:
    map_gen(ros::NodeHandle &nh_, int num_ob):nh(nh_)
    {
        num_obstacle = num_ob;
        Eigen::Vector3d mp_size(10.0, 10.0, 3.0);
        map_size = mp_size;
        obstacle_pub = nh.advertise<visualization_msgs::MarkerArray>("/homework2/obstacle",10);


    }

    inline void genObstacles();

    inline void visualizeObstacles();

    inline void getPotentialCostGrad(const Eigen::VectorXd xi, Eigen::VectorXd &g, double &f);

    inline void testGradPotential(const Eigen::VectorXd xi);

    inline bool smoothedL1(const double &x,
                          const double& mu,
                          double &f,
                          double &df);


};

inline void map_gen::genObstacles()
{
    std::vector<std::pair<Eigen::Vector3d, double>> ret_pair;
    srand(time(NULL));//设置随机数种子，使每次产生的随机序列不同
    for (int i = 0; i<num_obstacle; i++)
    {
        double random_x, random_y, random_z, random_r;
        

            random_x = map_size(0)*(rand() % (99 + 1) / (float)(99 + 1)) - map_size(0) /2.0;
            random_y = map_size(1)*(rand() % (99 + 1) / (float)(99 + 1)) - map_size(1) /2.0;
            //random_z = map_size(2)*(rand() % (99 + 1) / (float)(99 + 1)) - map_size(2) /2.0;
            random_z = 0.5;
            Eigen::Vector3d p_i (random_x, random_y, random_z);   
            random_r = 0.75 * (rand() % (99 + 1) / (float)(99 + 1));

            std::pair<Eigen::Vector3d, double> p1;
            p1 = std::make_pair(p_i, random_r);
            ret_pair.push_back(p1);          
         
    }

    obstacle_pair = ret_pair;

    
    return ;
}


inline void map_gen::visualizeObstacles()
{
    visualization_msgs::Marker obstacle_marker;
    visualization_msgs::MarkerArray ob_array;
    obstacle_marker.id = 0;
    obstacle_marker.type = visualization_msgs::Marker::SPHERE;
    
    obstacle_marker.header.frame_id = "map";
    obstacle_marker.pose.orientation.w = 1.00;
    obstacle_marker.action = visualization_msgs::Marker::ADD;
    obstacle_marker.ns = "traj";
    obstacle_marker.color.r = 1.00;
    obstacle_marker.color.g = 0.00;
    obstacle_marker.color.b = 0.00;
    obstacle_marker.color.a = 1.00;
    
    for (int i = 0; i< obstacle_pair.size(); i++)
    {
        obstacle_marker.id  += 1;
        Eigen::Vector3d pos = obstacle_pair[i].first;
        double r_i = obstacle_pair[i].second;
        obstacle_marker.header.stamp = ros::Time::now();
        obstacle_marker.pose.position.x = pos(0);
        obstacle_marker.pose.position.y = pos(1);
        obstacle_marker.pose.position.z = pos(2);
        obstacle_marker.scale.x = 2*r_i;
        obstacle_marker.scale.y = 2*r_i;
        obstacle_marker.scale.z= 2*r_i;
        ob_array.markers.push_back(obstacle_marker);

    }


    obstacle_pub.publish(ob_array);



}

inline void map_gen::getPotentialCostGrad(const Eigen::VectorXd xi, Eigen::VectorXd &g, double &f)
{
    int size_n = xi.size() /3;
    Eigen::VectorXd g_tmp  = Eigen::VectorXd::Zero(xi.size());
    int size_o = obstacle_pair.size();

    for (int i = 0; i<size_n;i++)
    {
        Eigen::Vector3d point_i(xi(3*i), xi(3*i+1), xi(3*i+2));
        //std::cout << "point_i in out_loop " << point_i.transpose() <<std::endl;
        Eigen::MatrixXd grad_xi_wrt_x = Eigen::MatrixXd::Zero(3,xi.size());
        grad_xi_wrt_x.block(0, 3*i, 3, 3) = Eigen::Matrix3d::Identity();

        for (int j = 0; j < size_o; j++)
        {
            //std::cout << "size_o: "  << size_o<< std::endl;
            Eigen::Vector3d obs_j = obstacle_pair[j].first;
            double r_j = obstacle_pair[j].second;

            //std::cout <<"obj111: " << obs_j.transpose() << std::endl;
            //std::cout << "r_j111: " <<  r_j << std::endl;


            double criterion_ij = r_j - (point_i - obs_j).norm();
            double smooth_factor = 0.01;
            double f_smooth = 0.0, g_smooth = 0.0;

            //if (smoothedL1(criterion_ij, smooth_factor, f_smooth, g_smooth ))
            if(criterion_ij >= 0)
            {
                if (print_flag)
                {
                    if (criterion_ij ==0)
                    {
                        std::cout <<"jjjjjjjjjjjjjjjjjjjjjj" << std::endl;
                    } 
                std::cout <<"obj: " << obs_j.transpose() << std::endl;
                std::cout << "r_j: " <<  r_j << std::endl;
                std::cout << "point_i: " << point_i.transpose() << std::endl;
                std::cout << "criterion_ij:" << criterion_ij << std::endl;

                }

                //f += 1000 * f_smooth;
                f += 1000 * criterion_ij;
                Eigen::VectorXd g_ij_wrt_x(xi.size());
                Eigen::Vector3d g_ij_wrt_xi;
                if((point_i - obs_j).norm() == 0.0)
                {
                    std::cout << "norm error================jjjjjjjjjjjjjjjjjjj===============" << std::endl;

                }
                g_ij_wrt_xi = -1000 *  0.5 * (1/ (point_i -obs_j).norm()) * (2 * point_i - 2* obs_j);
                g_ij_wrt_x = grad_xi_wrt_x.transpose() * g_ij_wrt_xi;

                //g_tmp +=  g_smooth * g_ij_wrt_x;
                g_tmp +=   g_ij_wrt_x;
                
            }
        }
        

    }
    g += g_tmp;
}


 inline void map_gen::testGradPotential(const Eigen::VectorXd xi)
 {

    print_flag = true;
    double delta = 0.00001;
    Eigen::VectorXd x_p_test = xi;
    Eigen::VectorXd x_m_test = xi;
    Eigen::VectorXd g_test(xi.size());
    for (int i = 0; i< xi.size(); i++)
    {
        x_p_test(i) += delta;
        double f_p_test=0.0;
        Eigen::VectorXd g_p_test(xi.size());
        getPotentialCostGrad(x_p_test,g_p_test,f_p_test);
        x_p_test(i) -= delta;
        x_m_test(i) -= delta;
        double f_m_test = 0.0;
        Eigen::VectorXd g_m_test(xi.size());
        getPotentialCostGrad(x_m_test,g_p_test,f_m_test);
        x_m_test(i) += delta;

        double g_test_i = (f_p_test - f_m_test) / (2* delta);

        g_test(i) = g_test_i;
     
    }

    Eigen::VectorXd g_analytic  = Eigen::VectorXd::Zero(xi.size());
    double f_analytic = 0;
    getPotentialCostGrad(xi, g_analytic, f_analytic);

    print_flag = false;

    std::cout << "grad_test============================================" << std::endl;
    std::cout << g_test.transpose() << std::endl;
    std::cout << "grad_analytical============================================" << std::endl;
    std::cout << g_analytic.transpose() << std::endl;



 }

inline bool map_gen::smoothedL1(const double &x,
                          const double& mu,
                          double &f,
                          double &df)
{
            if (x < 0.0)
            {
                return false;
            }
            else if (x > mu)
            {
                f = x - 0.5 * mu;
                df = 1.0;
                return true;
            }
            else
            {
                const double xdmu = x / mu;
                const double sqrxdmu = xdmu * xdmu;
                const double mumxd2 = mu - 0.5 * x;
                f = mumxd2 * sqrxdmu * xdmu;
                df = sqrxdmu * ((-0.5) * xdmu + 3.0 * mumxd2 / mu);
                return true;
            }

}
#endif



