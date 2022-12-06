#ifndef TOPP_H
#define TOPP_H

#include <Eigen/Eigen>
#include <iostream>
#include <algorithm>
#include <cmath>



class Topp
{
private:

//arc length parameters
 std::vector<double> s;
 std::vector<Eigen::Vector3d> qs;
 std::vector<Eigen::Vector3d> qds;
 std::vector<Eigen::Vector3d> qdds; 
 
 double b0 = 0.0, bl = 0.0;

 //objective function parameters
 Eigen::VectorXd c;  //cost function vector

//optimiazation variables
 Eigen::VectorXd xi; 
 Eigen::VectorXd grad_alm;
 Eigen::MatrixXd hess_alm;
 double f_xi;


 //constraint variables
 Eigen::MatrixXd G;  //G is the matrix of the equality constraints
 Eigen::VectorXd h;  //h is the vector of the equality constraints

 Eigen::MatrixXd A_linear_ineq_;
 Eigen::VectorXd b_linear_ineq_;

 std::vector<Eigen::MatrixXd> A_socp_ck_;
 std::vector<Eigen::Vector3d> b_socp_ck_;

 std::vector<Eigen::MatrixXd> A_socp_bk_;
 std::vector<Eigen::VectorXd> b_socp_bk_;


 //Augmented Lagrangian variables
    double rho;
    Eigen::VectorXd lambda;
    Eigen::VectorXd mu_ineq;
    std::vector<Eigen::Vector3d> mu_socp_ck;
    std::vector<Eigen::Vector3d> mu_socp_bk;

 //physical parameters
 double v_max = 3.0;
 double a_max = 3.0;


 //physical variables
 std::vector<Eigen::Vector3d> vel_vec_;
 std::vector<Eigen::Vector3d> acc_vec_;
 Eigen::VectorXd time_vec_;
 double duration_ = 0.0;
 double traj_start_time_ = 0.0;



public:
    inline void setData(const std::vector<double> &s_data,
                        const std::vector<Eigen::Vector3d> &qs_data,
                        const std::vector<Eigen::Vector3d> &qds_data,
                        const std::vector<Eigen::Vector3d> &qdds_data)
    {
        s.resize(s_data.size());
        qs.resize(qs_data.size());
        qds.resize(qds_data.size());
        qdds.resize(qdds_data.size());
        std::cout <<"444444444444444" << std::endl;
        for (int i = 0; i < s.size(); i++)
        {
            s[i] = s_data[i];
            qs[i] = qs_data[i];
            qds[i] = qds_data[i];
            qdds[i] = qdds_data[i];
        }
        // s = s_data;
        // this->qs = qs_data;
        // this->qds = qds_data;
        // this->qdds = qdds_data;
        std::cout <<"555555555555" << std::endl;
        xi = Eigen::VectorXd::Zero(s.size()*4);
        std::cout <<"666666666666" << std::endl;
    }
        
    

    inline void setConstraints();

    inline void evalCostAndGrad(const Eigen::VectorXd &x,
                                Eigen::VectorXd &grad,
                                double &f);

    inline void evalHession(const Eigen::VectorXd &x,
                       Eigen::MatrixXd &hess);

    

    inline Eigen::VectorXd projOnSecCone(const Eigen::VectorXd &v);
    inline Eigen::MatrixXd bSubDiffSOCP(const Eigen::VectorXd &x);

    inline int armijo_line_search(const Eigen::VectorXd &x,
                                const Eigen::VectorXd &grad_init,
                                const Eigen::VectorXd & dir,
                                double &f_init,
                                double &step
                                );

    inline int modified_damped_newton_optimiz();

    inline int solve();

    inline bool stopCriteria();

    inline void updateAugmentedLagrangianParameters();

    inline void testGrad(const Eigen::VectorXd &xi);
    inline void testHess(const Eigen::VectorXd &xi, Eigen::MatrixXd &hess);

    inline void setTimeVariables();
    
    inline void getPosVelAcc(Eigen::Vector3d & pos, Eigen::Vector3d& vel, Eigen::Vector3d& acc, double time);

    inline void setTrajStartTime(double time){traj_start_time_ = time;}

    inline void getDurationAndStartTime(double& duration, double& start_time){duration = duration_; start_time = traj_start_time_;}

   


    





};

    inline void Topp::setConstraints()
    {
        //set cost function vector
        std::cout <<"77777777777777" << std::endl;
        c = Eigen::VectorXd::Zero(s.size()*4);
        for (int i = 0; i<s.size()-1; i++)
        {
            c(s.size()*3+i) = 2*(s[i+1] - s[i]);
        }
        std::cout <<"8888888888888" << std::endl;

        //1. set G and g===============================================================
        G = Eigen::MatrixXd::Zero(s.size()+1,s.size()*4);
        h = Eigen::VectorXd::Zero(s.size()+1);
        for (int i = 0; i < s.size()-1; i++)
        {
            //b(k+1) -  bk = 2*[s(k+1) - s(k)] *ak;
            G(i,s.size()+i) = -1;
            G(i,s.size()+i+1) = 1;
            G(i,i) = -2*(s[i+1] - s[i]);
            //h(i) = 2* (s[i+1]-s[i]) * xi(i); 
            h(i) = 0.0;
        }
        G(s.size()-1,s.size()) = 1;
        h(s.size()-1) = b0;
        G(s.size(),2*s.size()-1) = 1;
        h(s.size()) = bl;

        std::cout <<"9999999" << std::endl;

        //2. set linear inequality constraints //===============================================================
        Eigen::MatrixXd A_linear_ineq_bk = Eigen::MatrixXd::Zero(s.size(),s.size()*4);
        Eigen::VectorXd b_linear_ineq_bk = Eigen::VectorXd::Zero(s.size());
        for (int i = 0; i < s.size(); i++)
        {
            A_linear_ineq_bk(i,s.size()+i) = -1;
        }
        
        std::cout <<"9999999111111111111111111111111" << std::endl;

        Eigen::MatrixXd A_linear_ineq_vk = Eigen::MatrixXd::Zero(s.size()*6,s.size()*4);
        Eigen::VectorXd b_linear_ineq_vk = Eigen::VectorXd::Zero(s.size()*6);
        for (int i = 0; i < s.size(); i++)
        {
            A_linear_ineq_vk(i*6,s.size()+i) = (qds[i](0))*(qds[i](0));
            b_linear_ineq_vk(i*6) = -v_max*v_max;
            A_linear_ineq_vk(i*6+1,s.size()+i) = -(qds[i](0))*(qds[i](0));
            b_linear_ineq_vk(i*6+1) = -v_max*v_max;
            A_linear_ineq_vk(i*6+2,s.size()+i) = (qds[i](1))*(qds[i](1));
            b_linear_ineq_vk(i*6+2) = -v_max*v_max;
            A_linear_ineq_vk(i*6+3,s.size()+i) = -(qds[i](1))*(qds[i](1));
            b_linear_ineq_vk(i*6+3) = -v_max*v_max;
            A_linear_ineq_vk(i*6+4,s.size()+i) = (qds[i](2))*(qds[i](2));
            b_linear_ineq_vk(i*6+4) = -v_max*v_max;
            A_linear_ineq_vk(i*6+5,s.size()+i) = -(qds[i](2))*(qds[i](2));
            b_linear_ineq_vk(i*6+5) = -v_max*v_max;
        }
        std::cout <<"99999991222222222222222222222" << std::endl;

        Eigen::MatrixXd A_linear_ineq_ak = Eigen::MatrixXd::Zero(s.size()*6,s.size()*4);
        Eigen::VectorXd b_linear_ineq_ak = Eigen::VectorXd::Zero(s.size()*6);
        for (int i = 0; i < s.size(); i++)
        {
            A_linear_ineq_ak(i*6,i) = (qds[i](0));
            A_linear_ineq_ak(i*6,i+s.size()) = (qdds[i](0));
            b_linear_ineq_ak(i*6) = -a_max;

            A_linear_ineq_ak(i*6+1,i) = -(qds[i](0));
            A_linear_ineq_ak(i*6+1,i+s.size()) = -(qdds[i](0));
            b_linear_ineq_ak(i*6+1) = -a_max;

            A_linear_ineq_ak(i*6+2,i) = (qds[i](1));
            A_linear_ineq_ak(i*6+2,i+s.size()) = (qdds[i](1));
            b_linear_ineq_ak(i*6+2) = -a_max;

            A_linear_ineq_ak(i*6+3,i) = -(qds[i](1));
            A_linear_ineq_ak(i*6+3,i+s.size()) = -(qdds[i](1));
            b_linear_ineq_ak(i*6+3) = -a_max;

            A_linear_ineq_ak(i*6+4,i) = (qds[i](2));
            A_linear_ineq_ak(i*6+4,i+s.size()) = (qdds[i](2));
            b_linear_ineq_ak(i*6+4) = -a_max;

            A_linear_ineq_ak(i*6+5,i) = -(qds[i](2));
            A_linear_ineq_ak(i*6+5,i+s.size()) = -(qdds[i](2));
            b_linear_ineq_ak(i*6+5) = -a_max;

        }
        std::cout <<"9999999133333333333333333333" << std::endl;

        Eigen::MatrixXd A_linear_ineq(s.size() + 6*s.size() + 6*s.size(),s.size()*4);
        A_linear_ineq << A_linear_ineq_bk,
                         A_linear_ineq_vk,
                         A_linear_ineq_ak;
        Eigen::VectorXd b_linear_ineq(s.size() + 6*s.size() + 6*s.size());
        b_linear_ineq << b_linear_ineq_bk,
                         b_linear_ineq_vk,
                         b_linear_ineq_ak;

        A_linear_ineq_ = A_linear_ineq;
        b_linear_ineq_ = b_linear_ineq;

        std::cout <<"100000" << std::endl;

        //3. set SOCP constraints //===============================================================
        std::vector<Eigen::MatrixXd> A_socp_ck;
        std::vector<Eigen::Vector3d> b_socp_ck;
        for (int i = 0; i < s.size()-1; i++)
        {
            Eigen::MatrixXd A_socp_ck_i = Eigen::MatrixXd::Zero(3,s.size()*4);
            Eigen::VectorXd b_socp_ck_i = Eigen::VectorXd::Zero(3);

            //cTx+d
            A_socp_ck_i(0,2*s.size()+i) = 1;
            A_socp_ck_i(0,2*s.size()+1 +i) = 1;
            A_socp_ck_i(0,3*s.size()+i) = 1;
            b_socp_ck_i(0) = 0;

            //Ax+b
            //A_socp_ck_i.row(1) = 0;
            b_socp_ck_i(1) = 2;

            A_socp_ck_i(2,2*s.size()+i) = 1;
            A_socp_ck_i(2,2*s.size()+1 +i) = 1;
            A_socp_ck_i(2,3*s.size()+i) = -1;
            b_socp_ck_i(2) = 0;

            A_socp_ck.push_back(A_socp_ck_i);
            b_socp_ck.push_back(b_socp_ck_i);
        }
        A_socp_ck_ = A_socp_ck;
        b_socp_ck_ = b_socp_ck;

        std::cout <<"1110000" << std::endl;
        std::vector<Eigen::MatrixXd> A_socp_bk;
        std::vector<Eigen::VectorXd> b_socp_bk;
        for (int i = 0; i < s.size(); i++)
        {
            Eigen::MatrixXd A_socp_bk_i = Eigen::MatrixXd::Zero(3,s.size()*4);
            Eigen::VectorXd b_socp_bk_i = Eigen::VectorXd::Zero(3);
            //cTx+d
            A_socp_bk_i(0,s.size()+i) = 1;
            b_socp_bk_i(0) = 1;

            //Ax+b
            A_socp_bk_i(1,2*s.size()+i) = 2;
            b_socp_bk_i(1) = 0;

            A_socp_bk_i(2,s.size()+i) = 1;
            b_socp_bk_i(2) = -1;
            A_socp_bk.push_back(A_socp_bk_i);
            b_socp_bk.push_back(b_socp_bk_i);
        }
        A_socp_bk_ = A_socp_bk;
        b_socp_bk_ = b_socp_bk;

        //Initialize the Augmented Lagrangian variables
        rho = 1.0;
        lambda = Eigen::VectorXd::Zero(h.size());
        mu_ineq = Eigen::VectorXd::Zero(b_linear_ineq_.size());
        mu_socp_ck = std::vector<Eigen::Vector3d>(b_socp_ck_.size(),Eigen::VectorXd::Zero(3));
        mu_socp_bk = std::vector<Eigen::Vector3d>(b_socp_bk_.size(),Eigen::VectorXd::Zero(3));



    }

    inline void Topp::evalCostAndGrad(const Eigen::VectorXd &x,
                                      Eigen::VectorXd &grad,
                                      double &f)
    {
       f = 0.0;
        Eigen::VectorXd grad_tmp = Eigen::VectorXd::Zero(x.size());

        //1. objective
        f+= c.dot(x) + 0.5 * rho * (G*x - h + lambda/rho).squaredNorm();
        grad_tmp += c + rho * G.transpose() * (G*x - h + lambda/rho);

        //2. inequality constraints
        Eigen::VectorXd linear_ineq_res = A_linear_ineq_ * x + b_linear_ineq_ + mu_ineq/rho;
        for (int i = 0; i < linear_ineq_res.size(); i++)
        {
            if (linear_ineq_res(i) > 0)
            {
                f += 0.5* rho * linear_ineq_res(i) * linear_ineq_res(i);
                grad_tmp += rho * A_linear_ineq_.row(i).transpose() * linear_ineq_res(i);
            }
        }

        //3. SOCP constraints-ck
        for (int i = 0; i < A_socp_ck_.size(); i++)
        {
            Eigen::VectorXd socp_ck_res = mu_socp_ck[i]/rho - A_socp_ck_[i] * x  - b_socp_ck_[i] ;
            Eigen::VectorXd v_ret = this->projOnSecCone(socp_ck_res);
            f += 0.5* rho * v_ret.squaredNorm();
            grad_tmp -= A_socp_ck_[i].transpose() * this->projOnSecCone(mu_socp_ck[i] - rho * (A_socp_ck_[i] * x + b_socp_ck_[i]));

        }

        //4. SOCP constraints-bk
        for (int i = 0; i < A_socp_bk_.size(); i++)
        {
            Eigen::VectorXd socp_bk_res = mu_socp_bk[i]/rho - A_socp_bk_[i] * x  - b_socp_bk_[i] ;
            Eigen::VectorXd v_ret = this->projOnSecCone(socp_bk_res);
            f += 0.5* rho * v_ret.squaredNorm();
            grad_tmp -= A_socp_bk_[i].transpose() * this->projOnSecCone(mu_socp_bk[i] - rho * (A_socp_bk_[i] * x + b_socp_bk_[i]));
        }

        grad = grad_tmp;
        

    }

    inline void Topp::testGrad(const Eigen::VectorXd &xi)
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
        evalCostAndGrad(x_p_test,g_p_test,f_p_test);
        //getCostandGrad(x_p_test,start, goal, g_p_test,f_p_test);
        x_p_test(i) -= delta;
        x_m_test(i) -= delta;
        double f_m_test = 0.0;
        Eigen::VectorXd g_m_test(xi.size());
        evalCostAndGrad(x_m_test,g_m_test,f_m_test);
        //getCostandGrad(x_m_test, start, goal, g_p_test,f_m_test);
        x_m_test(i) += delta;

        double g_test_i = (f_p_test - f_m_test) / (2* delta);

        g_test(i) = g_test_i;
     
    }

    Eigen::VectorXd g_analytic  = Eigen::VectorXd::Zero(xi.size());
    double f_analytic = 0;
    evalCostAndGrad(xi,g_analytic,f_analytic);

    //getCostandGrad(xi, start, goal, g_analytic, f_analytic);

    //print_flag = false;

    std::cout << "grad_test============================================" << std::endl;
    std::cout << g_test.transpose() << std::endl;
    std::cout << "grad_analytical============================================" << std::endl;
    std::cout << g_analytic.transpose() << std::endl;
    std::cout << "grad_diff============================================" << std::endl;
    std::cout << (g_test - g_analytic).transpose() << std::endl;
    std::cout <<"grad_diff_lifinity_norm============================================" << std::endl;
    //std::cout << (g_test - g_analytic).lpNorm<Eigen::Infinity>() << std::endl;
    std::cout << (g_test - g_analytic).lpNorm<Eigen::Infinity>() / xi.lpNorm<Eigen::Infinity>()<< std::endl;
}

inline void Topp::testHess(const Eigen::VectorXd &xi, Eigen::MatrixXd &hess)
{
    double delta = 0.0000001;
    Eigen::VectorXd x_p_test = xi;
    Eigen::VectorXd x_m_test = xi;
    Eigen::VectorXd g_test(xi.size());
    double f_tmp ;
    evalCostAndGrad(xi,g_test,f_tmp);
    Eigen::MatrixXd H_test = Eigen::MatrixXd::Zero(xi.size(),xi.size());
    for (int i = 0; i<g_test.size();i++)
    {
     
            x_p_test(i) += delta;
            x_m_test(i) -= delta;
            double f_p_test=0.0;
            double f_m_test = 0.0;
            Eigen::VectorXd g_p_test(xi.size());
            Eigen::VectorXd g_m_test(xi.size());
            evalCostAndGrad(x_p_test,g_p_test,f_p_test);
            evalCostAndGrad(x_m_test,g_m_test,f_m_test);
            x_p_test(i) -= delta;
            x_m_test(i) += delta;
            H_test.col(i) = (g_p_test - g_m_test) / (2* delta);

        
    }
    Eigen::MatrixXd H_analytic = Eigen::MatrixXd::Zero(xi.size(),xi.size());
    evalHession(xi,H_analytic);

    std::cout << "hess_diff============================================" << std::endl;
    std::cout << (H_test - H_analytic).lpNorm<Eigen::Infinity>() / xi.lpNorm<Eigen::Infinity>() << std::endl;
    std::cout << (H_test.transpose() - H_analytic).norm() << std::endl;

    std::cout << "hess_test - hess_analytic============================================" << std::endl;
    std::cout << (H_test - H_analytic).maxCoeff()<< std::endl;
    std::cout << (H_test - H_analytic).minCoeff()<< std::endl;
    std::cout << (H_test.transpose() - H_analytic).maxCoeff()<< std::endl;
    std::cout << (H_test.transpose() - H_analytic).minCoeff()<< std::endl;

    hess = H_test;







}

    inline void Topp::evalHession(const Eigen::VectorXd &x,
                       Eigen::MatrixXd &hess)
    {
        //1. hession of objective
        //std::cout <<"hession 1111111111" << std::endl;
        hess = Eigen::MatrixXd::Zero(x.size(),x.size());
        hess += rho * G.transpose() * G;

        //2. hession of inequality constraints
        Eigen::VectorXd  linear_ineq_res = A_linear_ineq_ * x + b_linear_ineq_ + mu_ineq/rho;
        Eigen::MatrixXd ambk_max_tmp  = Eigen::MatrixXd::Zero(linear_ineq_res.size(),linear_ineq_res.size());
        for (int i = 0; i< linear_ineq_res.size(); i++)
        {
            if (linear_ineq_res(i) >= 0.0)
            {
                hess += rho * A_linear_ineq_.row(i).transpose() * A_linear_ineq_.row(i);
                //ambk_max_tmp(i,i) = 1;
            }
        }
        //hess += rho * A_linear_ineq_.transpose() * ambk_max_tmp * A_linear_ineq_;
        //std::cout <<"hession 22222222222" << std::endl;
        //3. hession of SOCP constraints-ck
        for (int i = 0; i < A_socp_ck_.size(); i++)
        {
            Eigen::VectorXd socp_ck_res = mu_socp_ck[i] - rho * (A_socp_ck_[i] * x  + b_socp_ck_[i]);
            Eigen::MatrixXd abpk_tmp = this->bSubDiffSOCP(socp_ck_res);
            hess += rho * A_socp_ck_[i].transpose() * abpk_tmp * A_socp_ck_[i];
        }
        //std::cout <<"hession 33333333333" << std::endl;
        
        //4. hession of SOCP constraints-bk 
        for (int i = 0; i < A_socp_bk_.size(); i++)
        {
            Eigen::VectorXd socp_bk_res = mu_socp_bk[i] - rho * (A_socp_bk_[i] * x  + b_socp_bk_[i]);
            Eigen::MatrixXd abpk_tmp = this->bSubDiffSOCP(socp_bk_res);
            hess += rho * A_socp_bk_[i].transpose() * abpk_tmp * A_socp_bk_[i];
        }
        //std::cout <<"hession 44444444444" << std::endl;
      
    }



    inline Eigen::VectorXd Topp::projOnSecCone(const Eigen::VectorXd &v)
    {
        double v0 = v(0);
        Eigen::VectorXd v1 = v.tail(v.rows()-1);
        if (v0 <= -v1.norm())
        {
            return Eigen::VectorXd::Zero(v.rows());
        }
        if (std::fabs(v0) < v1.norm())
        {
            Eigen::VectorXd v1_tmp (v.rows());
            v1_tmp << v1.norm(), v1;
            Eigen::VectorXd v_tmp = (v0+v1.norm())/(2*v1.norm())*v1_tmp;
            return v_tmp;
        }
        else
        {
            Eigen::VectorXd v_tmp = v;
            return v_tmp;
        }

    }


    inline Eigen::MatrixXd Topp::bSubDiffSOCP(const Eigen::VectorXd &x)
    {
        //std::cout <<"bSubDiffSOCP: x_size" << x.size() << std::endl;
        double x1 = x(0);
        Eigen::VectorXd x2 = x.tail(x.rows()-1);
        //std::cout <<"bSubDiffSOCP: x2" << x2.size() << std::endl;
        if(x2.norm() <= x1)
        {
            return Eigen::MatrixXd::Identity(x.rows(),x.rows());
        }
        if (x2.norm() <= -x1)
        {
            return Eigen::MatrixXd::Zero(x.rows(),x.rows());
        }
        else
        {
            Eigen::MatrixXd A_tmp (x.rows(),x.rows());
            A_tmp(0,0)= 0.5;
            A_tmp.block(0,1,1,x.rows()-1) = x2.transpose()/(2*x2.norm());
            //std::cout << "5.44444444444444444" << std::endl;
            A_tmp.block(1,0,x.rows()-1,1) = x2/(2*x2.norm());
            //std::cout << "5.5555555" << std::endl;
            A_tmp.block(1,1,x.rows()-1,x.rows()-1) =(x1+x2.norm())/(2*x2.norm())*Eigen::MatrixXd::Identity(x.rows()-1,x.rows()-1)-x1*x2*x2.transpose()/(2.0 * std::pow(x2.norm(),3)); 
            //std::cout << "5.66666666666" << std::endl;
            return A_tmp;
        }

    }


    inline int Topp::armijo_line_search(const Eigen::VectorXd &x,
                                const Eigen::VectorXd &grad_init,
                                const Eigen::VectorXd & dir,
                                double &f_init,
                                double &step
                                )
    {
        int iter = 0;
        double dg_init = 0.0001 * grad_init.dot(dir);
        double f_new = f_init;
        Eigen::VectorXd grad_tmp;
        while (true)
        {
            iter++;
            Eigen::VectorXd x_new = x + step * dir;
            evalCostAndGrad(x_new,grad_tmp,f_new);
            //std::cout << "f_new: " << f_new << std::endl;
            //std::cout << "f_init: " << f_init << std::endl;
            //std::cout << "dg_init: " << dg_init << std::endl;
            if (f_new <= f_init + step * dg_init)
            {
                return iter;
            }
            step *= 0.5;
            if (step < 1e-8)
            {
                std::cout << "armijo line search failed: " << std::endl;
                return -1;
            }
            if (iter > 100)
            {
                std::cout << "armijo line search failed: maximum linesearch" << std::endl;
                return -2;
            }
        }

        return 0;
    }


    inline int Topp::modified_damped_newton_optimiz()
    { 
        int iter = 0;
        //std::cout <<"modified_damped_newton_optimiz: " << "11111111111111111111111111" << std::endl;
        while(true)
        {
            iter++;
            //std::cout << "444444444444" << std::endl;
            //std::cout <<"modified_damped_newton_optimiz: " << "2222222222222222222222" << std::endl;
            this->evalCostAndGrad(xi,grad_alm,f_xi);
            //std::cout <<"modified_damped_newton_optimiz: " << "333333333333333333333333333" << std::endl;
            //std::cout << "5555555555555" << std::endl;
            double inter_criterion = grad_alm.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>();
            std::cout <<"inter_criterion: " << inter_criterion << std::endl;
            if (inter_criterion < 1e-5)
            {
                return iter;
            }
            double epsilon = std::min(1.0,grad_alm.lpNorm<Eigen::Infinity>())/10;
            Eigen::MatrixXd hess(xi.rows(),xi.rows());
            //std::cout << "5.1111111111111111" << std::endl;
            //std::cout <<"modified_damped_newton_optimiz: " << "33333111111111111111111111111" << std::endl;
            evalHession(xi,hess);
            //testHess(xi, hess);
            //std::cout <<"modified_damped_newton_optimiz: " << "444444444444444444444" << std::endl;
            //std::cout << "5.2222222222222222" << std::endl;
            Eigen::MatrixXd hess_damped = hess + epsilon*Eigen::MatrixXd::Identity(hess.rows(),hess.cols());
            Eigen::LLT<Eigen::MatrixXd> llt(hess_damped);
            Eigen::VectorXd dir = llt.solve(-grad_alm);
            //std::cout <<"modified_damped_newton_optimiz: " << "555555555555555555555555" << std::endl;
        //std::cout << "6666666666666" << std::endl;
            double step = 1;
            int status = armijo_line_search(xi,grad_alm,dir,f_xi,step);
            if (status == -1)
            {
                return -1;
            }
            if (status == -2)
            {
                return -2;
            }
            xi = xi + step * dir;
            if (iter > 100)
            {
                //testGrad(xi);
                //testHess(xi,hess);
                std::cout <<" modified damped newton optimiz failed: maximum iterations" << std::endl;
                return -3;
                
            }

            std::cout <<"inner iter: " << iter << std::endl;
                    //<< "xi: " << xi.transpose() << std::endl;
                    
        }
    }


    inline int Topp::solve()
    {
        int iter = 0;
        while(true)
        {
            iter++;
            std::cout <<"topp outer==================================" << std::endl
            << "out iterations: " << iter << std::endl
            << "objective value: " << c.dot(xi) << std::endl;
           // << "xi: " << xi.transpose() << std::endl;
            //std::cout << "2222" << std::endl;
            if (this->stopCriteria())
            {

                //testGrad(xi);
                //Eigen::MatrixXd hess(xi.rows(),xi.rows());
                //testHess(xi,hess);
                std::cout <<"topp success========================" << std::endl
                        << "iterations: " << iter << std::endl
                        << "objective value: " << c.dot(xi) << std::endl;
                        //<< "xi: " << xi.transpose() << std::endl;
                return iter;
            }
            int status = modified_damped_newton_optimiz();

            //std::cout << "3333333333" << std::endl;
            if (status < 0 )
            {
                std::cout << "modified damped newton optimiz failed" << std::endl;
                return -1;
            }
            //if ((mu/rho - projOnSymCone(mu/rho - A*xi - b)).lpNorm<Eigen::Infinity>() < 1e-5 && grad_alm.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>() < 1e-6)

            if (iter > 100)
            {
                // testGrad(xi);
                // Eigen::MatrixXd hess(xi.rows(),xi.rows());
                // testHess(xi,hess);
                
                std::cout << ": maximum iterations" << std::endl;
                return -2;
            }
            //mu = projOnSymCone(mu/rho - A*xi - b);
            //rho = std::min((1+1)*rho,1e+3);
            this->updateAugmentedLagrangianParameters();

        
        }
        return 0;

    }

    inline bool Topp::stopCriteria()                             
    {
        double gp_norm = grad_alm.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>();

        double gxh_max = (G*xi - h).lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>();

        Eigen::VectorXd g_ieq = A_linear_ineq_ * xi + b_linear_ineq_;
        for (int i = 0; i< g_ieq.size(); i++)
        {
            g_ieq(i) = std::max(g_ieq(i),(-mu_ineq/rho)(i));
        }

        double g_ieq_max = g_ieq.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>();

        double socp_ck_max = 0.0;
        for (int i = 0; i< b_socp_ck_.size(); i++)
        {
            Eigen::VectorXd tmp_vec = (mu_socp_ck)[i]/rho - this->projOnSecCone(mu_socp_ck[i]/rho - A_socp_ck_[i]*xi - b_socp_ck_[i]);
            double tmp_max = tmp_vec.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>();

            socp_ck_max = std::max(socp_ck_max,tmp_max);
        }

        double socp_bk_max = 0.0;
        for (int i = 0; i< b_socp_bk_.size(); i++)
        {
            Eigen::VectorXd tmp_vec = (mu_socp_bk)[i]/rho - this->projOnSecCone(mu_socp_bk[i]/rho - A_socp_bk_[i]*xi - b_socp_bk_[i]);
            double tmp_max = tmp_vec.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>();

            socp_bk_max = std::max(socp_bk_max,tmp_max);
        }

        double kkt_max = std::max({gxh_max, g_ieq_max, socp_ck_max, socp_bk_max});

        // std::cout <<"outer criteria: " << std::endl
        // <<"gp_norm: " << gp_norm << std::endl
        // <<"kkkt_max: " << kkt_max << std::endl;

        if (gp_norm < 1e-6 && kkt_max < 1e-4)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    inline void Topp::updateAugmentedLagrangianParameters()
    {

        lambda = lambda + rho * (G* xi - h);

        Eigen::VectorXd mu_plus_rho_g = mu_ineq + rho * (A_linear_ineq_ * xi + b_linear_ineq_);
        for (int i = 0; i< b_linear_ineq_.size(); i++)
        {
    
            mu_ineq(i) = std::max(mu_plus_rho_g(i), 0.0);
        }

        for (int i = 0; i< b_socp_ck_.size(); i++)
        {
            mu_socp_ck[i] = this->projOnSecCone(mu_socp_ck[i] - rho*(A_socp_ck_[i]*xi + b_socp_ck_[i]));
        }

        for (int i = 0; i< b_socp_bk_.size(); i++)
        {
            mu_socp_bk[i] = this->projOnSecCone(mu_socp_bk[i] - rho*(A_socp_bk_[i]*xi + b_socp_bk_[i]));
        }

        rho = std::min((1+2)*rho,1e+3);


    }

    inline void Topp::setTimeVariables()
    {
        Eigen::VectorXd ak (s.size());
        Eigen::VectorXd bk (s.size());
        Eigen::VectorXd ck (s.size());
        Eigen::VectorXd dk (s.size());

        for (int i = 0; i< s.size(); i++)
        {
            ak(i) = xi(i);
            bk(i) = xi(s.size()+i);
            ck(i) = xi(2*s.size()+i);
            dk(i) = xi(3*s.size() + i);
        }

        Eigen::VectorXd time_vec  = Eigen::VectorXd::Zero(s.size());
        time_vec(0) = 0.0;
        double time_sum = 0.0;
        for (int i = 0; i< s.size()-1; i++)
        {
            time_sum += 2*(s[i+1] - s[i])*dk(i);
            time_vec(i+1) = time_sum;
            // std::cout << "time_vec: " << time_vec(i) << std::endl;
        }

        std::vector<Eigen::Vector3d> vel_vec(s.size());
        std::vector<Eigen::Vector3d> acc_vec(s.size());
        for (int i = 0; i< s.size(); i++)
        {
            vel_vec[i] = qds[i]*sqrt(bk(i));
            acc_vec[i] = qdds[i]*bk(i) + qds[i]*ak(i);

            // std::cout << "vel_vec: " << vel_vec[i].transpose() << std::endl;
            // std::cout << "acc_vec: " << acc_vec[i].transpose() << std::endl;


        }

        vel_vec_ = vel_vec;
        acc_vec_ = acc_vec;
        time_vec_ = time_vec;
        duration_ = time_sum;


        
    }

    inline void Topp::getPosVelAcc(Eigen::Vector3d & pos, Eigen::Vector3d& vel, Eigen::Vector3d& acc, double time)
    {
        if ( time >=0 && time < duration_)
        {
            int idx = 0;
        for (int i = 0; i< time_vec_.size(); i++)
        {
            if (time_vec_(i) > time)
            {
                idx = i;
                break;
            }
        }
        vel = vel_vec_[idx - 1] + (vel_vec_[idx] - vel_vec_[idx - 1])*(time - time_vec_(idx - 1))/(time_vec_(idx) - time_vec_(idx - 1));
        acc = acc_vec_[idx - 1] + (acc_vec_[idx] - acc_vec_[idx - 1])*(time - time_vec_(idx - 1))/(time_vec_(idx) - time_vec_(idx - 1));
        pos = qs[idx-1] + (qs[idx] - qs[idx-1])*(time - time_vec_(idx - 1))/(time_vec_(idx) - time_vec_(idx - 1));
   
        }
        else
        {
            vel = Eigen::Vector3d::Zero();
            acc = Eigen::Vector3d::Zero();

            
        }
 
    }












































#endif // TOPP_H