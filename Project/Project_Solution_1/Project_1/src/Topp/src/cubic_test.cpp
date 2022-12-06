#include <iostream>
#include <Eigen/Eigen>
#include "cubic_spline_banded_test.h"
//#include "cubic_spline.h"
#include <ros/ros.h>
#include "support/lbfgs.hpp"
#include "support/map_gen.h"
#include <geometry_msgs/PoseStamped.h>
#include "polymap.hpp"
#include "topp.h"
#include <std_msgs/Float64.h>
#include <visualization_msgs/Marker.h>

class cubic_test
{
private:
    Eigen::Vector3d start, goal;
    Eigen::VectorXd xi;
    ros::NodeHandle nh;
    Eigen::Vector3d mp_size;
    ros::Subscriber target_sub;
    std::vector<Eigen::Vector3d> startGoal;
    PolyMap polymap;
    cubic_spline cubic;
    map_gen  mp_gen;
    Topp topp_tester;

    //vellocity and acceleration publisher
    ros::Publisher vel_pub_0, vel_pub_1, vel_pub_2;
    ros::Publisher acc_pub_0, acc_pub_1, acc_pub_2;
    ros::Publisher spherePub;
    
public:
    
    cubic_test(Eigen::Vector3d p1, Eigen::Vector3d p2, ros::NodeHandle &nh_):cubic(nh_), mp_gen(nh_, 20), polymap(nh_, Eigen::Vector3d(10.0, 10.0, 10.0), 20)
    {
        start = p1;
        goal = p2;
        nh = nh_;
        target_sub = nh.subscribe("/move_base_simple/goal", 1, &cubic_test::targetCallback, this);
        vel_pub_0 = nh.advertise<std_msgs::Float64>("/vel_0", 1);
        vel_pub_1 = nh.advertise<std_msgs::Float64>("/vel_1", 1);
        vel_pub_2 = nh.advertise<std_msgs::Float64>("/vel_2", 1);
        acc_pub_0 = nh.advertise<std_msgs::Float64>("/acc_0", 1);
        acc_pub_1 = nh.advertise<std_msgs::Float64>("/acc_1", 1);
        acc_pub_2 = nh.advertise<std_msgs::Float64>("/acc_2", 1);
        spherePub = nh.advertise<visualization_msgs::Marker>("/sphere", 1);
        mp_gen.genObstacles();
        polymap.generate_obstacles();
        
    }

    int run(const int N)
    {
        int ret;
        if(startGoal.size() == 2)
        {
        start = startGoal[0];
        goal = startGoal[1];
        //mp_gen.genObstacles();
        double finalCost;
        Eigen::VectorXd x(3*N);

        Eigen::Vector3d delta_vec = (goal -start) /(N+1);

        /* Set the initial guess */
        for (int i = 0; i < N; i += 1)
        {

            x(3*i) = start(0)+delta_vec(0)*(i+1);
            x (3*i +1) = start(1)+delta_vec(1)*(i+1);
            x (3*i +2) = start(2)+delta_vec(2)*(i+1);
        }

        /* Set the minimization parameters */
        lbfgs::lbfgs_parameter_t params;
        params.g_epsilon = 1.0e-8;
        params.past = 3;
        params.delta = 1.0e-8;

        std::cout <<"init guess:: " << x.transpose() << std::endl;
        ros::Time t_start = ros::Time::now();

        /* Start minimization */
         ret = lbfgs::lbfgs_optimize(x,
                                        finalCost,
                                        cubiccost,
                                        nullptr,
                                        this,
                                        params);
        ros::Time t_end = ros::Time::now();

        double duration = (t_end - t_start).toSec();

        
        visualize();
        //cubic.testEnergyGrad(x, start, goal);
        //mp_gen.testGradPotential(x);

        //testFinalGrad(x, start, goal);

        /* Report the result. */
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "L-BFGS Optimization Returned: " << ret << std::endl
                  << "Minimized Cost: " << finalCost << std::endl
                  <<" Duration: " << duration << std::endl
                  << "Optimal Variables: " << std::endl
                  << x.transpose() << std::endl;

        // std::vector<Eigen::Vector3d> pos, vel, acc;
        // double dura_tmp = cubic.getDuration();
        // for(double s = 0.0; s <=dura_tmp; s += 0.01)
        // {
        // Eigen::Vector3d pos_tmp, vel_tmp, acc_tmp;
        // cubic.getPosVelAcc(s, pos_tmp, vel_tmp, acc_tmp);
        // pos.emplace_back(pos_tmp);
        // vel.emplace_back(vel_tmp);
        // acc.emplace_back(acc_tmp);
        // }
        // for(int i = 0; i < 300; i++)
        // {
        // std::cout <<"pos:"  << i << " "<<pos[i].transpose() <<std::endl;
        // std::cout <<"traj:" << i<< " " << cubic.traj_points[i].transpose() <<std::endl;
        // }
        cubic.setArcLenthPath(0.1);
        // for (int i = 0; i< cubic.arc_data_.s.size()-1; i++)
        // {
        //     std::cout  << std::setprecision(10) <<"index:" << i << "=========================================================== "<<std::endl  
        //     << "s: " << cubic.arc_data_.s[i] <<std::endl
        //     << "qs+1 - qs" << (cubic.arc_data_.qs[i+1] - cubic.arc_data_.qs[i]).transpose() <<std::endl
        //     << "qs" << cubic.arc_data_.qs[i].transpose() <<std::endl
        //     << "qds" << cubic.arc_data_.qds[i].transpose() <<std::endl
        //     << "qdds" << cubic.arc_data_.qdds[i].transpose() <<std::endl;

        // }

        
        // std::cout <<"111111" << std::endl;
        // std::cout <<"arc_data_.s.size(): " << cubic.arc_data_.s.size() << std::endl;
        // std::cout <<"arc_data_.qs.size(): " << cubic.arc_data_.qs.size() << std::endl;
        // std::cout <<"arc_data_.qds.size(): " << cubic.arc_data_.qds.size() << std::endl;
        // std::cout <<"arc_data_.qdds.size(): " << cubic.arc_data_.qdds.size() << std::endl;

        std::vector<double> s;
        std::vector<Eigen::Vector3d> qs, qds, qdds;

        for (int i = 0; i< cubic.arc_data_.s.size(); i += 30)
        {
            s.emplace_back(cubic.arc_data_.s[i]);
            qs.emplace_back(cubic.arc_data_.qs[i]);
            qds.emplace_back(cubic.arc_data_.qds[i]);
            qdds.emplace_back(cubic.arc_data_.qdds[i]);
        }

        // s.emplace_back(cubic.arc_data_.s.back());
        // qs.emplace_back(cubic.arc_data_.qs.back());
        // qds.emplace_back(cubic.arc_data_.qds.back());
        // qdds.emplace_back(cubic.arc_data_.qdds.back());

        // for (int i = 0; i< s.size(); i++)
        // {
        //     std::cout  << std::setprecision(10) <<"index:" << i << "=========================================================== "<<std::endl  
        //     << "s: " << s[i] <<std::endl
        //     << "qs" << qs[i].transpose() <<std::endl
        //     << "qds" << qds[i].transpose() <<std::endl
        //     << "qdds" << qdds[i].transpose() <<std::endl;
        // }
    
        //topp_tester.setData(cubic.arc_data_.s, cubic.arc_data_.qs, cubic.arc_data_.qds, cubic.arc_data_.qdds);
        topp_tester.setData(s, qs, qds, qdds);
        topp_tester.setConstraints();
        //std::cout <<"22222222222" << std::endl;
        topp_tester.solve();
        topp_tester.setTimeVariables();
        double time_start = ros::Time::now().toSec();
        //std::cout <<"333333" << std::endl;
        topp_tester.setTrajStartTime(time_start);

        double duration_topp , start_topp;
        topp_tester.getDurationAndStartTime(duration_topp, start_topp);
        std::cout <<"duration_topp: " << duration_topp << std::endl;



        }

        return ret;
    }

    void visualize()
    {
        //mp_gen.genObstacles();
        //mp_gen.visualizeObstacles();
        polymap.visualize_obstacles();
        cubic.visualizeSpline();
    }



    inline void targetCallback(const geometry_msgs::PoseStamped::ConstPtr &msg)
    {
        //std::cout <<"11" <<std::endl;
            if (startGoal.size() >= 2)
            {
                startGoal.clear();
            }
            const double zGoal = 0.5;
            const Eigen::Vector3d goal(msg->pose.position.x, msg->pose.position.y, zGoal);
            startGoal.emplace_back(goal);
            run(20);
            visualize();


            return;

    }

    inline void testFinalGrad(Eigen::VectorXd x, Eigen::Vector3d start, Eigen::Vector3d goal)
    {
        double delta = 0.00001;
        Eigen::VectorXd x_p, x_m;
        x_p = x;
        x_m = x;

        
        Eigen::VectorXd g_test(x.size());
        for  (int i = 0; i < x.size(); i++ )
        {
            x_p(i) += delta;
            Eigen::VectorXd gp (x.size());
            
            double f_p = 0.0;
            cubic.getCostandGrad(x_p, start, goal, gp, f_p);
            mp_gen.getPotentialCostGrad(x_p, gp, f_p);
            x_p(i) -= delta;

            x_m(i) -= delta;
            double f_m = 0.0;
            cubic.getCostandGrad(x_m, start, goal, gp, f_m);
            mp_gen.getPotentialCostGrad(x_m, gp, f_m);
            x_m(i) += delta;


            double g_i = (f_p - f_m ) / (2.0 * delta);

            g_test(i)  = g_i;
               
        }
        double f_analytic = 0.0;
        Eigen::VectorXd g_analytic = Eigen::VectorXd::Zero(x.size());
        cubic.getCostandGrad(x, start, goal, g_analytic, f_analytic);
        mp_gen.getPotentialCostGrad(x, g_analytic, f_analytic);

        std::cout << "final_g_test====================================" << std::endl;
        std::cout << g_test.transpose() << std::endl;
        std::cout <<"final_g_analytic=================================" << std::endl;
        std::cout << g_analytic.transpose() << std::endl;

    }

    inline void process()
    {
       //std::cout <<"process====================================" << std::endl;

        double duration , topp_traj_start_time;
        topp_tester.getDurationAndStartTime(duration, topp_traj_start_time);
        
        const double delta = ros::Time::now().toSec() - topp_traj_start_time;

        //std::cout <<"delta: " << delta << std::endl;
        //std::cout <<"duration: " << duration << std::endl;

        if(delta > 0.0 && delta < duration )
        {
            //std::cout <<"process22222222222222222222222222====================================" << std::endl;
            Eigen::Vector3d pos_tmp, vel_tmp, acc_tmp;
            topp_tester.getPosVelAcc(pos_tmp, vel_tmp, acc_tmp, delta);

            std_msgs::Float64 vel_msg_0, vel_msg_1, vel_msg_2, acc_msg_0, acc_msg_1, acc_msg_2;

            vel_msg_0.data = vel_tmp(0);
            vel_msg_1.data = vel_tmp(1);
            vel_msg_2.data = vel_tmp(2);
            acc_msg_0.data = acc_tmp(0);
            acc_msg_1.data = acc_tmp(1);
            acc_msg_2.data = acc_tmp(2);
            vel_pub_0.publish(vel_msg_0);
            vel_pub_1.publish(vel_msg_1);
            vel_pub_2.publish(vel_msg_2);
            acc_pub_0.publish(acc_msg_0);
            acc_pub_1.publish(acc_msg_1);
            acc_pub_2.publish(acc_msg_2);


            visualization_msgs::Marker sphereMarkers, sphereDeleter;

            sphereMarkers.id = 0;
            sphereMarkers.type = visualization_msgs::Marker::SPHERE_LIST;
            sphereMarkers.header.stamp = ros::Time::now();
            sphereMarkers.header.frame_id = "map";
            sphereMarkers.pose.orientation.w = 1.00;
            sphereMarkers.action = visualization_msgs::Marker::ADD;
            sphereMarkers.ns = "spheres";
            sphereMarkers.color.r = 1.00;
            sphereMarkers.color.g = 0.00;
            sphereMarkers.color.b = 0.00;
            sphereMarkers.color.a = 1.00;
            sphereMarkers.scale.x = 0.1 * 2.0;
            sphereMarkers.scale.y = 0.1 * 2.0;
            sphereMarkers.scale.z = 0.1 * 2.0;

            sphereDeleter = sphereMarkers;
            sphereDeleter.action = visualization_msgs::Marker::DELETE;

            geometry_msgs::Point point;
            point.x = pos_tmp(0);
            point.y = pos_tmp(1);
            point.z = pos_tmp(2);
            sphereMarkers.points.push_back(point);
            spherePub.publish(sphereDeleter);
            spherePub.publish(sphereMarkers);


        }


        
    }

    



    


private:
    static inline double cubiccost(void *instance,
                                  const Eigen::VectorXd &x,
                                 Eigen::VectorXd &g)
    {

        double f_ret = 0.0;
        cubic_test &obj = *(cubic_test *)instance;
        obj.cubic.getCostandGrad(x, obj.start, obj.goal, g, f_ret);

        //std::cout <<"energy cost::  " << f_ret << "============================="  <<std::endl;
        obj.polymap.getPotentialCostGrad(x, g, f_ret);
        
        //obj.mp_gen.getPotentialCostGrad(x, g, f_ret);
        //std::cout <<"potential cost:: " << f_ret << "===========================" << std::endl;

        //obj.mp_gen.testGradPotential(x);
        //std::cout << "g: " << " " << g.transpose() << std::endl;
        return f_ret;
    }

    
    static int monitorProgress(void *instance,
                               const Eigen::VectorXd &x,
                               const Eigen::VectorXd &g,
                               const double fx,
                               const double step,
                               const int k,
                               const int ls)
    {
        std::cout << std::setprecision(4)
                  << "================================" << std::endl
                  << "Iteration: " << k << std::endl
                  << "Function Value: " << fx << std::endl
                  << "Gradient Inf Norm: " << g.cwiseAbs().maxCoeff() << std::endl
                  << "Variables: " << std::endl
                  << x.transpose() << std::endl;
        return 0;
    }




};






int main(int argc, char **argv)
{
   
    ros::init(argc, argv, "cuibc_test_node");
    ros::NodeHandle nh_;
    //cubic_spline cubic(nh_);
    
    //Eigen::VectorXd xi(12);
    //xi << 2.0, 0.0, 1.0, 2.0, 4.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0;
    
    Eigen::Vector3d  start(0,0,1);
    Eigen::Vector3d goal(5,-4,1);


    
    //Eigen::Matrix<double, 12, -1> ret = cubic.getCoef(xi, start, goal);
    cubic_test cutest(start, goal, nh_);
    ros::Rate lr(1000);

    
    //cutest.run(4);
    //cutest.mp_gen.genObstacles();
    
    //ros::spinOnce();


    bool flag = true;
    int iter  = 0;
    
    while (ros::ok())
    {
        // if(flag)
        // {
        //     cutest.visualize();
        //     iter ++;
        //     if(iter > 10)
        //     {
        //         flag = false;
        //     }
            
        // }
        
        
        // std::vector<Eigen::Vector3d> trajpoints =  cubic.getTrajpoint(0.01, ret);
        // Eigen::VectorXd g(xi.size());
        // double f;
        // cubic.getCostandGrad(xi, start, goal, g, f);
        // cubic.visualizeSpline(trajpoints);
        //cutest.visualize();
        cutest.visualize();
        cutest.process();
        ros::spinOnce();
        lr.sleep();
    }

    return 0;
}