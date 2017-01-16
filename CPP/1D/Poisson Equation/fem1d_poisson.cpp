/* Solving 1D Poisson's equation with finite element method

	-u''(x) = f(x)	on	0 < x < 1
			with	u(0) = u(1) = 0

	Input:  m, number of nodes
	Output: u, FE solution of nodal point

*/ 

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/LU>
#include <cmath>

Eigen::MatrixXd FormStiffness(int);
Eigen::VectorXd FormRHS(int);
double Integrate_hat1(double, double);
double Integrate_hat2(double, double);
double psi1(double, double, double);
double psi2(double, double, double);
double f(double);


int main(int argc, char* argv[]) {
	
	argc = 1;
	int m = atoi(argv[1]);	// Number of nodes entered in prompt
	
	Eigen::MatrixXd K(m, m);
	K = FormStiffness(m);	// Constructing stiffness matrix
	
	Eigen::VectorXd b(m);
	b = FormRHS(m);		// Force vector

	Eigen::VectorXd u(m);
	u = K.lu().solve(b);	// Solving the system of equations Au = b
	std::cout << u << std::endl;		
	
	return 0;
	
}

Eigen::VectorXd FormRHS(int m) {
	
	double h = 1.0/(m - 1);
	double  x[m];

	Eigen::VectorXd b(m);	

	for (int i = 0; i < m; i++) {
		x[i] = i * h;
	}
	
	b(0) 	 = 0.0; // Account for BC start
	b(m - 1) = 0.0; // Account for BC end

	for (int j = 1; j < m - 1; j++) {
		b(j) = Integrate_hat1(x[j - 1], x[j]) + Integrate_hat2(x[j], x[j + 1]);
	}

	return b;
}

Eigen::MatrixXd FormStiffness(int m) {
	
	double h = 1.0/(m - 1);

	Eigen::MatrixXd K(m, m );
	K.setZero(m, m);

	for (int i = 1; i < m - 1; i++)
		K(i, i) = 2.0;
	for (int i = 1; i < m - 2; i++) {
		K(i, i + 1) = -1.0;
		K(i + 1, i) = -1.0;
	}
	
	K(0, 0) 	= 1.0; // Account for BC start
	K(m - 1, m - 1) = 1.0; // Account for BC end
	
	K *= 1.0/h;

	return K;
}

double Integrate_hat1(double a, double b) {
	
	double c = (a + b)/2.0;                                                                                                    
	return ((b - a)/6.0)*(f(a)*psi1(a, a, b) + 4*f(c)*psi1(c, a, b) + f(b)*psi1(b, a, b));	
}

double Integrate_hat2(double a, double b) {
	
	double c = (a + b)/2.0;
	return ((b - a)/6.0)*(f(a)*psi2(a, a, b) + 4*f(c)*psi2(c, a, b) + f(b)*psi2(b, a, b));
}

double psi1(double x, double a, double b) {
	
	double h = b - a;
	return (x - a)/h;
}

double psi2(double x, double a, double b) {
	
	double h = b - a;
	return (b - x)/h;
}


double f(double x) {

	return x*x - 2*x;
}



				



