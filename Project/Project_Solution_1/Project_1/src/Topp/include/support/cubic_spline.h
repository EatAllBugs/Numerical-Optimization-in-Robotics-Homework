#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include<Eigen/Eigen>
#include<visualization_msgs/Marker.h>
#include <geometry_msgs/Point.h>
#include<std_msgs/ColorRGBA.h>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <ros/ros.h>
#include "support/banded_system.h"

class cubic_spline
{
private:
    /* data */
    ros::NodeHandle nh;
    ros::Publisher spline_pub;
    std::vector<Eigen::Vector3d> traj_points;
    BandedSystem A;

    
public:


    cubic_spline(ros::NodeHandle &nh_):nh(nh_)
    {
        
        spline_pub = nh.advertise<visualization_msgs::Marker>("/homework2/spline", 10);
    }

    inline Eigen::Matrix<double, 12, -1> getCoef(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal);

    inline void visualizeSpline();
    
    inline std::vector<Eigen::Vector3d> getTrajpoint(double delta_t,Eigen::Matrix<double, 12, -1> coef);

    inline void getCostandGrad(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal, Eigen::VectorXd &g, double &f);

    inline void testEnergyGrad(const Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal);


    // inline double getFandGrad(Eigen::VectorXd xi, )

   
};

inline Eigen::Matrix<double, 12, -1>  cubic_spline::getCoef(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal)
{

    if (xi.size() % 3 != 0)
    {
        std::cout << "error: incorrect decesion variable " <<std::endl;

        return {};
    }
    int size_n = xi.size() /3 ;

    Eigen::Matrix3Xd x_all;
    x_all.resize(3,size_n +2);
    
    x_all.leftCols(1) = start;

    for(int i = 0; i< size_n; i++)
    {
        Eigen::Vector3d x_i(xi(3*i), xi(3*i+1), xi(3*i+2));
        x_all.col(i+1) = x_i;
    }
    x_all.rightCols(1) = goal;
    //std::cout << x_all << std::endl;

    
    
  
    Eigen::VectorXd b_right(3*size_n);
    Eigen::MatrixXd A_left = Eigen::MatrixXd::Zero(3*size_n, 3*size_n);

    Eigen::Matrix3d f_ele = 4*Eigen::Matrix3d::Identity();
    Eigen::Matrix3d unity_ele = Eigen::Matrix3d::Identity();



    for (int i = 0; i< size_n ; i++)
    {
        b_right.segment(3*i,3) = 3*(x_all.col(i+2)- x_all.col(i));
        if(i==0)
        {
            A_left.block(0,0,3,3)  = f_ele;
            A_left.block(0,3,3,3)  = unity_ele;  
            continue;
        }
        if(i == (size_n-1))
        {
            A_left.block(3*i, 3*size_n -6, 3, 3) = unity_ele;
            A_left.block(3*i, 3*size_n -3, 3, 3) = f_ele;
            continue;
        }
        A_left.block(3*i,3*(i-1),3,3) = unity_ele;
        A_left.block(3*i,3*(i),3,3) = f_ele;
        A_left.block(3*i,3*(i+1),3,3) = unity_ele;
   
    }
    // std::cout <<" bright: " <<  b_right.transpose() << std::endl;
    // std::cout <<"1111111111111111" <<std::endl;
    Eigen::MatrixXd A_inverse = A_left.inverse();

    //std::cout <<"A:==================================================" << std::endl; 

    //std::cout << A_left<< std::endl;

    //std::cout <<"A_inverse:========================================================" << std::endl;
    //std::cout <<A_inverse<< std::endl;

    Eigen::VectorXd coef_vec_Di = A_left.colPivHouseholderQr().solve(b_right);
    Eigen::VectorXd coef_vec_D_final;
    coef_vec_D_final.resize(coef_vec_Di.size()+6);
    coef_vec_D_final << 0.0, 0.0, 0.0, coef_vec_Di, 0.0, 0.0, 0.0;
    Eigen::Matrix <double, 12, -1> matrix_ret_coef;
    matrix_ret_coef.resize(12, size_n +1);

    for (int i = 0; i<size_n+1; i++)
    {
        Eigen::Vector3d Di = coef_vec_D_final.segment(3*i,3);
        Eigen::Vector3d Di1 = coef_vec_D_final.segment(3*(i+1),3);
        Eigen::Vector3d ai, bi, ci, di;
        ai = x_all.col(i);
        bi = Di;
        ci = 3*(x_all.col(i+1)- x_all.col(i)) - 2* Di - Di1;
        di = 2*(x_all.col(i)- x_all.col(i+1))+ Di +Di1;
        Eigen::VectorXd coef_i(12);
        coef_i << ai,bi,ci,di;
        matrix_ret_coef.col(i) = coef_i;

        // std::cout <<"xi: " << ai.transpose() <<std::endl;
        // std::cout <<"Di "  << Di.transpose() << std::endl;
        // std::cout <<"Di1 " << Di1.transpose() << std::endl;
        
    }

    return matrix_ret_coef;



    
}

inline std::vector<Eigen::Vector3d> cubic_spline::getTrajpoint(double delta_t, Eigen::Matrix<double, 12, -1> coef)
{
    std::vector<Eigen::Vector3d> points_ret;
    int size = coef.cols();
    for (int i = 0; i<size; i++)
    {
        Eigen::Vector3d ai, bi, ci, di;
        ai = coef.col(i).head(3);
        bi = coef.col(i).segment(3,3);
        ci = coef.col(i).segment(6,3);
        di = coef.col(i).tail(3);

        for (double t = 0.0; t<=1.01; t+=delta_t)
        {
            double t1 = t;
            double t2 = t1 * t1;
            double t3 = t2*t1;

            Eigen::Vector3d pos_i = ai + bi * t1 + ci * t2 + di *t3;
            // std::cout << " t " << t << std::endl; 
            // std::cout <<" iter point: " << pos_i.transpose()<< std::endl;
            points_ret.emplace_back(pos_i);
        }

    }

    return points_ret;
}

inline void cubic_spline::visualizeSpline()
{



    visualization_msgs::Marker  traj_marker;


        traj_marker.id = 0;
        traj_marker.type = visualization_msgs::Marker::LINE_LIST;
        traj_marker.header.stamp = ros::Time::now();
        traj_marker.header.frame_id = "map";
        traj_marker.pose.orientation.w = 1.00;
        traj_marker.action = visualization_msgs::Marker::ADD;
        traj_marker.ns = "traj";
        traj_marker.color.r = 0.00;
        traj_marker.color.g = 0.00;
        traj_marker.color.b = 1.00;
        traj_marker.color.a = 1.00;
        traj_marker.scale.x = 0.10;


        if(traj_points.size() == 0)
        {

            //std::cout << "[Cubic_Spline:]" << " visualizer: " << " empty trajpoints" << std::endl;
            return;
        }


        for (int i = 0; i< traj_points.size()-1; i++)
        {
            geometry_msgs::Point point;
            point.x = traj_points[i](0);
            point.y = traj_points[i](1);
            point.z = traj_points[i](2);

            traj_marker.points.push_back(point);
            
            point.x = traj_points[i+1](0);
            point.y = traj_points[i+1](1);
            point.z = traj_points[i+1](2);

            traj_marker.points.push_back(point);

            

        }

        spline_pub.publish(traj_marker);



}

inline void cubic_spline::getCostandGrad(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal, Eigen::VectorXd &g, double &f)
{



    Eigen::VectorXd g_tmp = Eigen::VectorXd::Zero(g.size());
    if (xi.size() % 3 != 0)
    {
        std::cout << "error: incorrect decesion variable " <<std::endl;

        return ;
    }
    int size_n = xi.size() /3 ;

    Eigen::Matrix3Xd x_all;
    x_all.resize(3,size_n +2);
    
    x_all.leftCols(1) = start;

    for(int i = 0; i< size_n; i++)
    {
        Eigen::Vector3d x_i(xi(3*i), xi(3*i+1), xi(3*i+2));
        x_all.col(i+1) = x_i;
    }
    x_all.rightCols(1) = goal;
    //std::cout << x_all << std::endl;

    
    
  
    Eigen::VectorXd b_right(3*size_n);
    Eigen::MatrixXd A_left = Eigen::MatrixXd::Zero(3*size_n, 3*size_n);

    Eigen::Matrix3d f_ele = 4*Eigen::Matrix3d::Identity();
    Eigen::Matrix3d unity_ele = Eigen::Matrix3d::Identity();



    for (int i = 0; i< size_n ; i++)
    {
        b_right.segment(3*i,3) = 3*(x_all.col(i+2)- x_all.col(i));
        if(i==0)
        {
            A_left.block(0,0,3,3)  = f_ele;
            A_left.block(0,3,3,3)  = unity_ele;  
            continue;
        }
        if(i == (size_n-1))
        {
            A_left.block(3*i, 3*size_n -6, 3, 3) = unity_ele;
            A_left.block(3*i, 3*size_n -3, 3, 3) = f_ele;
            continue;
        }
        A_left.block(3*i,3*(i-1),3,3) = unity_ele;
        A_left.block(3*i,3*(i),3,3) = f_ele;
        A_left.block(3*i,3*(i+1),3,3) = unity_ele;
   
    }
    // std::cout <<" bright: " <<  b_right.transpose() << std::endl;
    // std::cout <<"1111111111111111" <<std::endl;
    Eigen::MatrixXd A_inverse = A_left.inverse();

    //std::cout <<"A:==================================================" << std::endl; 

    //std::cout << A_left<< std::endl;

    //std::cout <<"A_inverse:========================================================" << std::endl;
    //std::cout <<A_inverse<< std::endl;

    Eigen::MatrixXd grad_bright_wrt_x  = Eigen::MatrixXd::Zero(3* size_n, 3* size_n);
    Eigen::Matrix3d ele_3 = 3.0*Eigen::Matrix3d::Identity();

    grad_bright_wrt_x.block(0,3,3,3) = ele_3;
    grad_bright_wrt_x.block(3*size_n-3, 3*size_n-6,3,3) = -ele_3;

    for (int i = 1; i<size_n-1; i++)
    {
        grad_bright_wrt_x.block(3*i, 3*(i-1), 3, 3) = -ele_3;
        grad_bright_wrt_x.block(3*i, 3*(i+1), 3, 3) = ele_3;    
    }

    //std::cout << "grad_bright_wrt_x=============================================" << std::endl;
    //std::cout << grad_bright_wrt_x << std::endl;


    Eigen::MatrixXd D0 = Eigen::MatrixXd::Zero(3,3*size_n);
    Eigen::MatrixXd DN = Eigen::MatrixXd::Zero(3,3*size_n);

    const int size_tmp = size_n;
    std::vector<Eigen::MatrixXd> grad_D_wrt_x;
    grad_D_wrt_x.push_back(D0);
    for (int i = 0; i<size_n; i++)
    {
        Eigen::MatrixXd grad_Di_wrt_x = A_inverse.block(3*i, 0, 3, 3*size_n)*grad_bright_wrt_x;
        grad_D_wrt_x.push_back(grad_Di_wrt_x);
        
    }
    grad_D_wrt_x.push_back(DN);



    Eigen::VectorXd coef_vec_Di = A_left.colPivHouseholderQr().solve(b_right);
    Eigen::VectorXd coef_vec_D_final;
    coef_vec_D_final.resize(coef_vec_Di.size()+6);
    coef_vec_D_final << 0.0, 0.0, 0.0, coef_vec_Di, 0.0, 0.0, 0.0;
    Eigen::Matrix <double, 12, -1> matrix_ret_coef;
    matrix_ret_coef.resize(12, size_n +1);

    


    for (int i = 0; i<size_n+1; i++)
    {
        Eigen::Vector3d Di = coef_vec_D_final.segment(3*i,3);
        Eigen::Vector3d Di1 = coef_vec_D_final.segment(3*(i+1),3);
        Eigen::Vector3d ai, bi, ci, di;
        ai = x_all.col(i);
        bi = Di;
        ci = 3*(x_all.col(i+1)- x_all.col(i)) - 2* Di - Di1;
        di = 2*(x_all.col(i)- x_all.col(i+1))+ Di +Di1;
        Eigen::VectorXd coef_i(12);
        coef_i << ai,bi,ci,di;
        matrix_ret_coef.col(i) = coef_i;

        double f_i  = 12*di.dot(di) + 12* ci.dot(di) + 4*ci.dot(ci);
        f += f_i;


        Eigen::MatrixXd grad_xi_wrt_x = Eigen::MatrixXd::Zero(3, 3*size_n);
        Eigen::MatrixXd grad_xi1_wrt_x = Eigen::MatrixXd::Zero(3, 3*size_n);

        if(i==0)
        {
            grad_xi1_wrt_x.block(0,3*i,3,3) = Eigen::Matrix3d::Identity();
            
            
        }
        else if( i == size_n)
        {
            grad_xi_wrt_x.block(0, 3*size_n-3, 3, 3) = Eigen::Matrix3d::Identity();
        }
        else
        {
            grad_xi_wrt_x.block(0, 3*(i-1), 3, 3) = Eigen::Matrix3d::Identity();
            grad_xi1_wrt_x.block(0, 3*i, 3, 3) = Eigen::Matrix3d::Identity();

        }

        Eigen::MatrixXd grad_Di_wrt_x = grad_D_wrt_x[i];
        Eigen::MatrixXd grad_Di1_wrt_x = grad_D_wrt_x[i+1];

        Eigen::MatrixXd grad_ci_wrt_x = 3* grad_xi1_wrt_x - 3* grad_xi_wrt_x - 2* grad_Di_wrt_x - grad_Di1_wrt_x;
        Eigen::MatrixXd grad_di_wrt_x = 2 * grad_xi_wrt_x - 2* grad_xi1_wrt_x + grad_Di_wrt_x + grad_Di1_wrt_x;

        Eigen::VectorXd gi = 24* grad_di_wrt_x.transpose()*di + 12 * grad_ci_wrt_x.transpose() * di + 
                            12* grad_di_wrt_x.transpose() * ci + 8 * grad_ci_wrt_x.transpose()*ci;

        g_tmp += gi;

        // std::cout <<"xi: " << ai.transpose() <<std::endl;
        // std::cout <<"Di "  << Di.transpose() << std::endl;
        // std::cout <<"Di1 " << Di1.transpose() << std::endl;      
    }

    std::vector<Eigen::Vector3d>  traj_points_tmp = getTrajpoint(0.01, matrix_ret_coef);
    traj_points = traj_points_tmp;
    //visualizeSpline(traj_points_tmp);

    g = g_tmp;



    




}

inline void cubic_spline::testEnergyGrad(const Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal)
{
    //print_flag = true;
    double delta = 0.0000001;
    Eigen::VectorXd x_p_test = xi;
    Eigen::VectorXd x_m_test = xi;
    Eigen::VectorXd g_test(xi.size());
    for (int i = 0; i< xi.size(); i++)
    {
        x_p_test(i) += delta;
        double f_p_test=0.0;
        Eigen::VectorXd g_p_test(xi.size());
        getCostandGrad(x_p_test,start, goal, g_p_test,f_p_test);
        x_p_test(i) -= delta;
        x_m_test(i) -= delta;
        double f_m_test = 0.0;
        Eigen::VectorXd g_m_test(xi.size());
        getCostandGrad(x_m_test, start, goal, g_p_test,f_m_test);
        x_m_test(i) += delta;

        double g_test_i = (f_p_test - f_m_test) / (2* delta);

        g_test(i) = g_test_i;
     
    }

    Eigen::VectorXd g_analytic  = Eigen::VectorXd::Zero(xi.size());
    double f_analytic = 0;
    getCostandGrad(xi, start, goal, g_analytic, f_analytic);

    //print_flag = false;

    std::cout << "grad_test_energy============================================" << std::endl;
    std::cout << g_test.transpose() << std::endl;
    std::cout << "grad_analytical_energy============================================" << std::endl;
    std::cout << g_analytic.transpose() << std::endl;
}


































#endif