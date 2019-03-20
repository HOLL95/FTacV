#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/pybind11.h>
#include </users/henney/Documents/Oxford/C++_libraries/pybind11/include/pybind11/stl.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>
#include <math.h>
#include <iostream>
#include <exception>
#include <vector>
namespace py = pybind11;
using namespace std;
struct e_surface_fun {
    double E,dE;
    double Cdl,CdlE,CdlE2,CdlE3;
    double E0;
    double Ru,R;
    double k0;
    double alpha;
    double In0,u1n0;
    double dt;
    double gamma;

    double exp11,exp12;
    double dexp11,dexp12;
    double u1n1;
    double du1n1;
    double Cdlp;

    e_surface_fun (
                    const double E,
                    const double dE,
                    const double Cdl,
                    const double CdlE,
                    const double CdlE2,
                    const double CdlE3,
                    const double E0,
                    const double Ru,
                    const double R,
                    const double k0,
                    const double alpha,
                    const double In0,
                    const double u1n0,
                    const double dt,
                    const double gamma

                    ) :
        E(E),dE(dE),Cdl(Cdl),CdlE(CdlE),CdlE2(CdlE2),CdlE3(CdlE3),E0(E0),Ru(Ru),R(R),k0(k0),alpha(alpha),In0(In0),u1n0(u1n0),dt(dt),gamma(gamma) { }

    boost::math::tuple<double,double> operator()(const double In1) {
        update_temporaries(In1);
        return boost::math::make_tuple(residual(In1),residual_gradient(In1));
    }

    double residual(const double In1) const {
        return Cdlp*(dt*dE-Ru*(In1-In0)) + dt*R*(E-Ru*In1) - dt*In1 + gamma*(u1n1-u1n0);
        //return Cdlp*(dt*dE) - dt*In1 + (u1n1-u1n0) + Ru*E*dt;
    }
    double residual_gradient(const double In1) const {
        return -Cdlp*Ru - dt*R*Ru - dt + gamma*du1n1;
        //return -Cdlp*Ru - dt + du1n1;
    }

    void update_temporaries(const double In1) {
        const double Ereduced = E - Ru*In1;
        //const double Ereduced = E;
        const double Ereduced2 = pow(Ereduced,2);
        const double Ereduced3 = Ereduced*Ereduced2;
        const double expval1 = Ereduced - E0;
        exp11 = std::exp((1.0-alpha)*expval1);
        exp12 = std::exp(-alpha*expval1);

        dexp11 = -Ru*(1.0-alpha)*exp11;
        dexp12 = Ru*alpha*exp11;

        const double u1n1_top = dt*k0*exp11 + u1n0;
        const double du1n1_top = dt*k0*dexp11;
        const double denom = (dt*k0*exp11 + dt*k0*exp12 + 1);
        const double ddenom = dt*k0*(dexp11 + dexp12);
        const double tmp = 1.0/denom;
        const double tmp2 = pow(tmp,2);
        u1n1 = u1n1_top*tmp;
        du1n1 = -(u1n1_top*ddenom + du1n1_top*denom)*tmp2;

        Cdlp = Cdl*(1.0 + CdlE*Ereduced + CdlE2*Ereduced2 + CdlE3*Ereduced3);
        //Cdlp = Cdl*(1.0 + CdlE*Edc+ CdlE2*pow(Edc,2)+ CdlE3*pow(Edc,3));
    }
};
double et(double E_start, double E_reverse, double omega, double phase, double delta_E,double t){
	double E_t=(E_start+delta_E)+delta_E*(std::sin((omega*t)+phase));

	return E_t;
}

double dEdt(double E_start, double E_reverse, double omega, double phase, double delta_E, double t){
  double dedt=(delta_E*omega*std::cos(omega*t+phase));
	return dedt;
}
double c_et(double E_start, double E_reverse, double omega, double phase, double v, double delta_E, double t){
	double E_dc;
	double E_t;
  double tr=((E_reverse-E_start)/v);
	if (t<tr){
		E_dc=E_start+(v*t);
	}else {
		E_dc=E_reverse-(v*(t-tr));
	}

	 E_t=E_dc+(delta_E*(sin((omega*t)+phase)));

	return E_t;
}

double c_dEdt( double E_start, double E_reverse, double omega, double phase, double v, double delta_E, double t){ //
	double E_dc;
  double tr=((E_reverse-E_start)/v);
	if (t < tr){
		 E_dc=v;
	}else {
		 E_dc=-v;
	}
  double dedt= E_dc+(delta_E*omega*cos(omega*t+phase));

	return dedt;
}
py::object e_surface(double Cdl, double CdlE, double CdlE2, double CdlE3,double omega,double v,double alpha ,double E_start, double E_reverse, \
	double delta_E, double Ru, int Nt,  std::vector<double> times, double gamma, double E0, double k0, double phase, double pi, int num_points) {
  //set up temporal mesh
  vector<double> Itot;
  const double dt = (1.0/Nt)*2*pi/omega;
  if (times.size()==0) {
      const double Tmax = 2*(E_reverse-E_start)/v;
      const int Nt = Tmax/dt;
      std::cout << "\tNt= "<<Nt<<std::endl;
      Itot.resize(Nt,0);
      times.resize(Nt);
      for (int i=0; i<Nt; i++) {
          times[i] = i*dt;
      }
  } else {
#ifndef NDEBUG
        //std::cout << "\thave "<<times.size()<<" samples from "<<times[0]<<" to "<<times[times.size()-1]<<std::endl;
#endif
        Itot.resize(times.size(),0);
    }

    double R=0;
    double Itot0,Itot1;
    double u1n0;
    double t1 = 0.0;
    u1n0 = 1.0;

    //const double E = c_et( E_start,  E_reverse,  omega,  phase,  v,  delta_E, t1);
    //const double dE = c_dEdt(E_start,  E_reverse,  omega,  phase,  v,  delta_E, t1+0.5*dt);
    const double E = et(E_start, E_reverse, omega, phase,delta_E , t1+dt);
    const double dE =dEdt(E_start, E_reverse, omega, phase,delta_E , t1+0.5*dt);
    const double Cdlp = Cdl*(1.0 + CdlE*E + CdlE2*pow(E,2)+ CdlE3*pow(E,2));
    const double Itot_bound = 100000;//std::max(10*Cdlp*delta_E*omega/Nt,1.0);
    //std::cout << "Itot_bound = "<<Itot_bound<<std::endl;
    const int digits_accuracy = std::numeric_limits<double>::digits;
    const boost::uintmax_t maxit = 50;
    Itot0 = Cdlp*dE;
    Itot1 = Itot0;
    for (int n_out = 0; n_out < times.size(); n_out++) {
        while (t1 < times[n_out]) {
            Itot0 = Itot1;
            const double E = et(E_start, E_reverse, omega, phase,delta_E ,t1+dt);
            const double dE =dEdt(E_start, E_reverse, omega, phase, delta_E ,t1+0.5*dt);
            //const double E = c_et( E_start,  E_reverse,  omega,  phase,  v,  delta_E, t1);
            //const double dE = c_dEdt(E_start,  E_reverse,  omega,  phase,  v,  delta_E, t1+0.5*dt);
            boost::uintmax_t it = maxit;
            e_surface_fun bc(E,dE,Cdl,CdlE,CdlE2,CdlE3,E0,Ru,R,k0,alpha,Itot0,u1n0,dt,gamma);
            Itot1 = boost::math::tools::newton_raphson_iterate(bc, Itot0,Itot0-Itot_bound,Itot0+Itot_bound, digits_accuracy);
            //if (it == maxit)
            //{
            //  throw std::runtime_error("non-linear solve for Itot[n+1] failed, max number of iterations reached");
            //}
            bc.update_temporaries(Itot1);
            u1n0 = bc.u1n1;
            t1 += dt;
        }
        Itot[n_out] = (Itot1-Itot0)*(times[n_out]-t1+dt)/dt + Itot0;

        //if (n_out < 5) {
        //    std::cout << "at n_out = "<<n_out<<" Itot = "<<Itot[n_out]<<Itot1<<std::endl;
        //}
    }
  return py::cast(Itot);
}
PYBIND11_MODULE(isolver_noramp, m) {
	m.def("e_surface", &e_surface, "solve for I_tot with dispersion");
  m.def("c_et", &c_et, "Classical potential input");
  m.def("c_dEdt", &c_dEdt, "Classical potential input derivative");
  m.def("et", &et, "Ramp-free potential input");
  m.def("dEdt", &c_et, "Ramp-free potential derivative");


}
