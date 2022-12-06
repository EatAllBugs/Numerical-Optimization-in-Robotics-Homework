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


struct arc_data
{
    std::vector<double> s;
    std::vector<Eigen::Vector3d> qs;
    std::vector<Eigen::Vector3d> qds;
    std::vector<Eigen::Vector3d> qdds;

};


class cubic_spline
{
private:
    /* data */
    ros::NodeHandle nh;
    ros::Publisher spline_pub;
    
    Eigen::Matrix <double, 12, -1> cubic_coef;


    BandedSystem A;
    Eigen::MatrixXd b;
    int N; // inter_point_size;
    Eigen::Vector3d start_pt, goal_pt;



    
public:
    std::vector<Eigen::Vector3d> traj_points;
    arc_data arc_data_;



    cubic_spline(ros::NodeHandle &nh_):nh(nh_)
    {
        
        spline_pub = nh.advertise<visualization_msgs::Marker>("/homework2/spline", 10);
    }

    inline Eigen::Matrix<double, 12, -1> getCoef(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal);

    inline void visualizeSpline();
    
    inline std::vector<Eigen::Vector3d> getTrajpoint(double delta_t,Eigen::Matrix<double, 12, -1> coef);

    inline void getCostandGrad(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal, Eigen::VectorXd &g, double &f);

    inline void testEnergyGrad(const Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal);

    //inline void getPosVelAcc(double s, Eigen::Vector3d &pos, Eigen::Vector3d &vel, Eigen::Vector3d &acc);

    inline double getDuration(){return double (N+1);}

    inline void setArcLenthPath(double delta_t);


    // inline double getFandGrad(Eigen::VectorXd xi, )

    inline void setConditions(const int & point_n )
    {


        N= point_n;
        A.create(3*N, 3, 3);
        b.resize(3*N, 1);
        // start_pt = start;
        // goal_pt = goals;

        return ;

    }

    inline void setParameters(const Eigen::VectorXd &xi, const Eigen::Vector3d &start, Eigen::Vector3d const &goal)
    {

        if (xi.size() % 3 != 0)
        {
            std::cout << "error: incorrect decesion variable " <<std::endl;

            return ;
        }

        int size_n = xi.size() /3 ;

        setConditions(size_n);



        Eigen::Matrix3Xd x_all;
        x_all.resize(3,size_n +2);

        x_all.leftCols(1) = start;

        for(int i = 0; i< size_n; i++)
        {
        Eigen::Vector3d x_i(xi(3*i), xi(3*i+1), xi(3*i+2));
        x_all.col(i+1) = x_i;
        }
        x_all.rightCols(1) = goal;

        
        
        
        // Eigen::Vector3Xd x_all;
        // x_all.resize(3, inPs.size()+2);
        // x_all.leftCols(1) = start_pt;
        // x_all.block(1,0,3,inPs.size()) = inPs;
        // x_all.rightCols(1) = goal_pt;

        A.reset();
        b.setZero();
        

        for(int i = 0; i< N; i++)
        {
            b(3*i,0) = 3*(x_all.col(i+2)-x_all.col(i))(0);
            b(3*i+1,0) = 3*(x_all.col(i+2)-x_all.col(i))(1);
            b(3*i+2,0) = 3*(x_all.col(i+2)-x_all.col(i))(2);
            if(i==0)
            {
                A(0,0) = 4;
                A(1,1) = 4;
                A(2,2) = 4;
                A(0,3) = 1;
                A(1,4) = 1;
                A(2,5) = 1;
                continue;
            }
            if (i == N-1)
            {
                A(3*i, 3*N-6) = 1;
                A(3*i, 3*N-3) = 4;
                A(3*i+1, 3*N-5) = 1;
                A(3*i+1, 3*N-2) = 4;
                A(3*i+2, 3*N-4) = 1;
                A(3*i+2, 3*N-1) = 4;
                continue;
            }

            A(3*i,3*(i-1)) = 1;
            A(3*i,3*i) = 4;
            A(3*i,3*(i+1)) = 1;
            A(3*i+1, 3*(i-1)+1) = 1;
            A(3*i+1, 3*(i)+1) = 4;
            A(3*i+1, 3*(i+1)+1) = 1;
            A(3*i+2, 3*(i-1)+2) = 1;
            A(3*i+2, 3*(i)+2) = 4;
            A(3*i+2, 3*(i+1)+2) = 1;
        }

        A.factorizeLU();
        A.solve(b);

        
    }

   
};

inline void cubic_spline::setArcLenthPath(double delta_t)
{
   double arc_lenth = 0;
   std::vector<double> arc_lenth_vec;
   std::vector<Eigen::Vector3d> qs_vec;
   std::vector<Eigen::Vector3d> qds_vec;
   std::vector<Eigen::Vector3d> qdds_vec;
   


   //initialize the q(s) at 0;
   arc_lenth_vec.push_back(0.0);
   Eigen::Vector3d qs0 = cubic_coef.col(0).head(3);
   qs_vec.push_back(qs0);

   for(double t = delta_t; t<=double(N+1); t+=delta_t)
   {
         Eigen::Vector3d pos_tmp, vel_tmp, acc_tmp;
         
        int j= int (t/1);
        double tj = t-j;
        if (tj < 0.9999*delta_t ) //tj == 0
        {
            j = j-1;
            tj = 1.0;

        }
        double tj2 = tj *tj;
        double tj3 = tj2 * tj;
        
        Eigen::Vector3d ai, bi, ci, di;
        ai = cubic_coef.col(j).head(3);
        bi = cubic_coef.col(j).segment(3,3);
        ci = cubic_coef.col(j).segment(6,3);
        di = cubic_coef.col(j).tail(3);

        pos_tmp = ai + bi*tj + ci*tj2 + di*tj3;
        

        // integration of the arc_length between tj-delta_t and tj 
        for (double t_int = tj - delta_t; t_int< 0.99999*tj; t_int += delta_t/10.0)
        {
            
            if (t_int < 0)
            {
                //std::cout << "error: t_int < 0" << std::endl;
                t_int = 0.0;
            }
            Eigen::Vector3d vel_tmp_low = bi + 2*ci*t_int + 3* di* t_int * t_int;
            double f_lower = sqrt(vel_tmp_low.squaredNorm()) ;

            double t_up = t_int + delta_t/10;
            Eigen::Vector3d vel_tmp_up = bi + 2* ci*t_up + 3* di* t_up * t_up;
            double f_up = sqrt(vel_tmp_up.squaredNorm()) ;

            arc_lenth += (f_lower + f_up)*delta_t/20;
            arc_lenth_vec.push_back(arc_lenth);

            Eigen::Vector3d qs_tmp = ai +bi*t_up +ci*t_up*t_up +di*t_up*t_up*t_up;
            qs_vec.push_back(qs_tmp);

           // std::cout << std::setprecision(6)<<"arc_sample :    t_up :" << t_up <<  "t_int:" << " " <<t_int << " tj: " << tj << " t:  " <<  t << std::endl;
            //std::cout <<"arc_sample :    qs_tmp :" << qs_tmp.transpose() << std::endl;
            //std::cout << "0.1 < 0.1 ?" << (0.1 < 0.1) << std::endl;

            
        }
   }
   qds_vec.resize(qs_vec.size());
   qdds_vec.resize(qs_vec.size());

   for (int i = 0; i<qs_vec.size(); i++)
   {
       if (i==0)
       {
        qds_vec[i] = (qs_vec[i+1]-qs_vec[i])/(arc_lenth_vec[i+1]-arc_lenth_vec[i]);
       }
         else if (i==qs_vec.size()-1)
         {
          qds_vec[i] = (qs_vec[i]-qs_vec[i-1])/(arc_lenth_vec[i]-arc_lenth_vec[i-1]);
         }
         else
         {
          qds_vec[i] = (qs_vec[i+1]-qs_vec[i-1])/(arc_lenth_vec[i+1]-arc_lenth_vec[i-1]);
         }
   }

   for (int i = 0; i<qs_vec.size(); i++)
   {
       if (i==0)
       {
        qdds_vec[i] = (qds_vec[i+1]-qds_vec[i])/(arc_lenth_vec[i+1]-arc_lenth_vec[i]);
       }
         else if (i==qs_vec.size()-1)
         {
          qdds_vec[i] = (qds_vec[i]-qds_vec[i-1])/(arc_lenth_vec[i]-arc_lenth_vec[i-1]);
         }
         else
         {
          qdds_vec[i] = (qds_vec[i+1]-qds_vec[i-1])/(arc_lenth_vec[i+1]-arc_lenth_vec[i-1]);
         }
   }

   arc_data_.s = arc_lenth_vec;
   arc_data_.qs = qs_vec;
   arc_data_.qds = qds_vec;
   arc_data_.qdds = qdds_vec;


}

inline Eigen::Matrix<double, 12, -1>  cubic_spline::getCoef(Eigen::VectorXd xi, Eigen::Vector3d start, Eigen::Vector3d goal)
{

    setParameters(xi, start, goal);

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
    // //std::cout << x_all << std::endl;

    
    
  
    // Eigen::VectorXd b_right(3*size_n);
    // Eigen::MatrixXd A_left = Eigen::MatrixXd::Zero(3*size_n, 3*size_n);

    // Eigen::Matrix3d f_ele = 4*Eigen::Matrix3d::Identity();
    // Eigen::Matrix3d unity_ele = Eigen::Matrix3d::Identity();



    // for (int i = 0; i< size_n ; i++)
    // {
    //     b_right.segment(3*i,3) = 3*(x_all.col(i+2)- x_all.col(i));
    //     if(i==0)
    //     {
    //         A_left.block(0,0,3,3)  = f_ele;
    //         A_left.block(0,3,3,3)  = unity_ele;  
    //         continue;
    //     }
    //     if(i == (size_n-1))
    //     {
    //         A_left.block(3*i, 3*size_n -6, 3, 3) = unity_ele;
    //         A_left.block(3*i, 3*size_n -3, 3, 3) = f_ele;
    //         continue;
    //     }
    //     A_left.block(3*i,3*(i-1),3,3) = unity_ele;
    //     A_left.block(3*i,3*(i),3,3) = f_ele;
    //     A_left.block(3*i,3*(i+1),3,3) = unity_ele;
   
    // }
    // // std::cout <<" bright: " <<  b_right.transpose() << std::endl;
    // // std::cout <<"1111111111111111" <<std::endl;
    // Eigen::MatrixXd A_inverse = A_left.inverse();

    // //std::cout <<"A:==================================================" << std::endl; 

    // //std::cout << A_left<< std::endl;

    // //std::cout <<"A_inverse:========================================================" << std::endl;
    // //std::cout <<A_inverse<< std::endl;

    // Eigen::VectorXd coef_vec_Di = A_left.colPivHouseholderQr().solve(b_right);
    // Eigen::VectorXd coef_vec_D_final;
    // coef_vec_D_final.resize(coef_vec_Di.size()+6);
    // coef_vec_D_final << 0.0, 0.0, 0.0, coef_vec_Di, 0.0, 0.0, 0.0;
    Eigen::Matrix <double, 12, -1> matrix_ret_coef;
    matrix_ret_coef.resize(12, N +1);


    Eigen::MatrixXd D_all;

    //Eigen::MatrixX3d D_all;
    Eigen::Vector3d D0 = Eigen::Vector3d::Zero();
    D_all.resize(b.rows()+6,1);
    D_all.topRows(3) = D0;
    D_all.block(3,0,b.rows(),1) = b;
    D_all.bottomRows(3) = D0;
    Eigen::VectorXd D_all_vec = D_all;

    for (int i = 0; i<N+1; i++)
    {
        
        Eigen::Vector3d Di = D_all_vec.segment(3*i,3);
        Eigen::Vector3d Di1 = D_all_vec.segment(3*(i+1),3);
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

        for (double t = 0.0; t<1.00; t+=delta_t)
        {
            //std::cout <<"t: " << t << std::endl;
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
    // //std::cout << x_all << std::endl;

    
    
  
     Eigen::VectorXd b_right(3*size_n);
    // Eigen::MatrixXd A_left = Eigen::MatrixXd::Zero(3*size_n, 3*size_n);

    // Eigen::Matrix3d f_ele = 4*Eigen::Matrix3d::Identity();
    // Eigen::Matrix3d unity_ele = Eigen::Matrix3d::Identity();



    // for (int i = 0; i< size_n ; i++)
    // {
    //     b_right.segment(3*i,3) = 3*(x_all.col(i+2)- x_all.col(i));
    //     if(i==0)
    //     {
    //         A_left.block(0,0,3,3)  = f_ele;
    //         A_left.block(0,3,3,3)  = unity_ele;  
    //         continue;
    //     }
    //     if(i == (size_n-1))
    //     {
    //         A_left.block(3*i, 3*size_n -6, 3, 3) = unity_ele;
    //         A_left.block(3*i, 3*size_n -3, 3, 3) = f_ele;
    //         continue;
    //     }
    //     A_left.block(3*i,3*(i-1),3,3) = unity_ele;
    //     A_left.block(3*i,3*(i),3,3) = f_ele;
    //     A_left.block(3*i,3*(i+1),3,3) = unity_ele;
   
    // }
    // // std::cout <<" bright: " <<  b_right.transpose() << std::endl;
    // // std::cout <<"1111111111111111" <<std::endl;
    // Eigen::MatrixXd A_inverse = A_left.inverse();

    //std::cout <<"A:==================================================" << std::endl; 

    //std::cout << A_left<< std::endl;

    //std::cout <<"A_inverse:========================================================" << std::endl;
    //std::cout <<A_inverse<< std::endl;

    setParameters(xi, start, goal);



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
    Eigen::MatrixXd grad_D_wrt_x_banded = grad_bright_wrt_x;
    A.solve(grad_D_wrt_x_banded);

    Eigen::MatrixXd D0 = Eigen::MatrixXd::Zero(3,3*size_n);
    Eigen::MatrixXd DN = Eigen::MatrixXd::Zero(3,3*size_n);

    const int size_tmp = size_n;
    std::vector<Eigen::MatrixXd> grad_D_wrt_x;
    grad_D_wrt_x.push_back(D0);
    for (int i = 0; i<size_n; i++)
    {
        Eigen::MatrixXd grad_Di_wrt_x = grad_D_wrt_x_banded.block(3*i, 0, 3, 3*size_n);
        grad_D_wrt_x.push_back(grad_Di_wrt_x);
        
    }
    grad_D_wrt_x.push_back(DN);



    Eigen::VectorXd coef_vec_Di = b;
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
    cubic_coef = matrix_ret_coef;
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