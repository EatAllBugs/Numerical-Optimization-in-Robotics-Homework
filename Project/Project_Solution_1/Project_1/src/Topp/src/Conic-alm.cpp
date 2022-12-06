# include <Eigen/Eigen>
# include <iostream>


class ConicAlm
{
private:

//Optimization variables
Eigen::VectorXd xi;
Eigen::VectorXd grad_alm;
double f_xi;
double stp;

//Objective function variables
Eigen::VectorXd c;

//SOCP CONSTRAINT VARIABLES
Eigen::MatrixXd A;
Eigen::VectorXd b;


//Augmented Lagrangian variables
Eigen::VectorXd mu;
double rho;








public:

void inline setCondition(Eigen::MatrixXd A1, Eigen::VectorXd b1,Eigen::VectorXd c1)
{
    this->A = A1;
    this->b = b1;
    this->c = c1;
    mu = Eigen::VectorXd::Zero(A.rows());
    rho = 1;
    xi = Eigen::VectorXd::Zero(A.cols());
}


void inline evalCostAndGrad(const Eigen::VectorXd &x,
                      Eigen::VectorXd &grad,
                      double &f)
{
    Eigen::VectorXd v = projOnSymCone(mu/rho - A*x - b);
    f  = c.dot(x) + 0.5*rho*v.squaredNorm();
    grad = c - A.transpose()*projOnSymCone(mu - rho*(A*x+b));

}

void inline evalHession(const Eigen::VectorXd &x,
                       Eigen::MatrixXd &hess)
{
    //std::cout << "5.33333333333333" << std::endl;
    Eigen::MatrixXd hess_temp = bSubDiffSOCP(mu - rho*(A*x+b));
    //std::cout << "5.jjjjjjjjjj" << std::endl;
    hess = rho*A.transpose()*hess_temp*A;
    //std::cout << "5.zzzzzzzzzzz" << std::endl;
}


int inline armijo_line_search(const Eigen::VectorXd &x,
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


Eigen::VectorXd inline projOnSymCone(const Eigen::VectorXd &v)
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

Eigen::MatrixXd inline bSubDiffSOCP(const Eigen::VectorXd &x)
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
        A_tmp.block(1,1,x.rows()-1,x.rows()-1) =(x1+x2.norm())/(2*x2.norm())*Eigen::MatrixXd::Identity(x.rows()-1,x.rows()-1)-x1*x2*x2.transpose()/std::pow(x2.norm(),3); 
        //std::cout << "5.66666666666" << std::endl;
        return A_tmp;
    }

}



int inline modified_damped_newton_optimiz()
{ 
    int iter = 0;
    while(true)
    {
        iter++;
        //std::cout << "444444444444" << std::endl;
        evalCostAndGrad(xi,grad_alm,f_xi);
        //std::cout << "5555555555555" << std::endl;
        if (grad_alm.norm()/xi.norm() < 1e-6)
        {
            return iter;
        }
        double epsilon = std::min(1.0,grad_alm.lpNorm<Eigen::Infinity>())/10;
        Eigen::MatrixXd hess(xi.rows(),xi.rows());
        //std::cout << "5.1111111111111111" << std::endl;
        evalHession(xi,hess);
        //std::cout << "5.2222222222222222" << std::endl;
        Eigen::MatrixXd hess_damped = hess + epsilon*Eigen::MatrixXd::Identity(hess.rows(),hess.cols());
        Eigen::LLT<Eigen::MatrixXd> llt(hess_damped);
        Eigen::VectorXd dir = llt.solve(-grad_alm);
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
            std::cout <<" modified damped newton optimiz failed: maximum iterations" << std::endl;
            return -3;
            
        }

        std::cout <<"inner iter: " << iter << std::endl
                  << "xi: " << xi.transpose() << std::endl;
                  
    }
}

int inline solve()
{
    int iter = 0;
    while(true)
    {
        iter++;
        std::cout <<"ConicALM outer==================================" << std::endl
        << "out iterations: " << iter << std::endl
        << "objective value: " << c.dot(xi) << std::endl
        << "xi: " << xi.transpose() << std::endl;
        //std::cout << "2222" << std::endl;
        int status = modified_damped_newton_optimiz();
        //std::cout << "3333333333" << std::endl;
        if (status < 0 )
        {
            std::cout << "modified damped newton optimiz failed" << std::endl;
            return -1;
        }
        if ((mu/rho - projOnSymCone(mu/rho - A*xi - b)).lpNorm<Eigen::Infinity>() < 1e-5 && grad_alm.lpNorm<Eigen::Infinity>()/xi.lpNorm<Eigen::Infinity>() < 1e-6)
        {
            std::cout <<"ConicALM success========================" << std::endl
                      << "iterations: " << iter << std::endl
                      << "objective value: " << c.dot(xi) << std::endl
                      << "xi: " << xi.transpose() << std::endl;
           return iter;
        }
        if (iter > 100)
        {
            std::cout << "modified damped newton optimiz failed: maximum iterations" << std::endl;
            return -2;
        }
        mu = projOnSymCone(mu/rho - A*xi - b);
        rho = std::min((1+1)*rho,1e+3);

       
    }
    return 0;

}



};

int main()
{

    Eigen::MatrixXd A(8,7);
    Eigen::VectorXd b(8);
    Eigen::VectorXd c(7);
    A << 1,0,0,0,0,0,0,
         7,0,0,0,0,0,0,
         0,6,0,0,0,0,0,
         0,0,5,0,0,0,0,
         0,0,0,4,0,0,0,
         0,0,0,0,3,0,0,
         0,0,0,0,0,2,0,
         0,0,0,0,0,0,1;
    b << 1,1,3,5,7,9,11,13;
    c << 1,2,3,4,5,6,7;
    ConicAlm conic_alm_tester;
    conic_alm_tester.setCondition(A,b,c);
    //std::cout << "1" << std::endl;
    int status = conic_alm_tester.solve();
    std::cout << "status: " << status << std::endl;
    return 0;

}