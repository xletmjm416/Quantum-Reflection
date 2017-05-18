#include "QMSystem.h"

int main() {
	//data input
	double step_t;
	double step_x;
	std::cout << "Input step in space: ";
	std::cin >> step_x;
	std::cout << "Input step in time: ";
	std::cin >> step_t;
	int N_space = (int)(1/step_x); //number of space points
	int N_time = (int)(1/step_t); //number of time points
	std::cout << "Total steps in time: " << N_time << std::endl;
	std::cout << "Total steps in space: " << N_space << std::endl;
	std::cout << "dx/dt: " << step_x/step_t << std::endl;
	
	//system init
	VectorC space = utils::linspace(step_x);
	VectorC psi = space.unaryExpr(&gauss_state);
	VectorR pot = VectorR::Zero(N_space);
	
	QMSystem system(psi, pot, step_x, step_t);
	std::cout << system.get_state().cwiseAbs2() << std::endl;
	//output init
	std::ofstream output;
	output.open("out.csv");
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", ";"); //copied from reference
	
	VectorR prob_distr = psi.cwiseAbs2();
	output << prob_distr.transpose() << std::endl;
	for(int t=0; t<N_time*10;t++) {
		psi = system.cranknicolson();
		prob_distr = psi.cwiseAbs2();
		std::cout << prob_distr.transpose() << std::endl;
		output << prob_distr.transpose().format(CommaInitFmt) << std::endl;
	}
	std::cout << "finished";
	output.close();
	char c;
	std::cin >> c;
	return 0;
}