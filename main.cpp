#include "main.h"

#define SYSTEM_SIZE 10
int main() {

	double h_bar=1, mass=1;
	//data input
	std::cout << "-------- Properties of space and time --------" << std::endl << std::endl;
	double step_t, step_x;
	std::cout << "Input step in space: ";
	std::cin >> step_x;
	std::cout << "Input step in time: ";
	std::cin >> step_t;
	int N_space = (int)(SYSTEM_SIZE/step_x); //number of space points
	int N_time = (int)(1/step_t); //number of time points
	std::cout << "System spatial length: " << SYSTEM_SIZE << std::endl << std::endl;
	std::cout << "Total steps in time: " << N_time << std::endl;
	std::cout << "Total steps in space: " << N_space << std::endl;
	std::cout << "Ratio step_x/step_t: " << step_x/step_t << std::endl << std::endl;
	std::cout << std::endl;
	
	VectorC space = utils::linspace(step_x, SYSTEM_SIZE);
	VectorC psi(N_space);

	std::cout << "-------- Initial parameters of Gaussian wavepacket --------" << std::endl << std::endl;
	double momentum, initial_pos, spread;
	std::cout << "Initial momentum: ";
	std::cin >> momentum;
	std::cout << "Initial position: ";
	std::cin >> initial_pos;
	std::cout << "Initial spread: ";
	std::cin >> spread;
	std::cout << std::endl;
	for(int i=0; i<N_space; i++) {
		dcplx val = utils::gauss_state(space(i), momentum, initial_pos, spread);
		psi(i) = val;
	}
	std::cout << "Velocity: " << h_bar*momentum/mass << std::endl << std::endl;
	
	std::cout << "-------- Potential function --------" << std::endl << std::endl;
	std::cout << "(not yet implemented)" << std::endl << std::endl;
	VectorR pot = VectorR::Zero(N_space);

	
	//system init
	QMSystem system(psi, pot, step_x, step_t, h_bar, mass, 1);
	
	//std::cout << system.get_state().cwiseAbs2() << std::endl;

	//output init
	std::ofstream output_x, output_p;
	output_x.open("out-space.csv");
	output_p.open("out-momentum.csv");
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ";", "", "", "", "", ""); //copied from reference
	
	output_x << psi.transpose().format(CommaInitFmt) << std::endl;
	VectorC psi_momentum = FFT::FFT(psi);
	output_p << psi_momentum.transpose().format(CommaInitFmt) << std::endl;
	for(int t=0; t<10*N_time;t++) {
		psi = system.cranknicolson();
		output_x << psi.transpose().format(CommaInitFmt) << std::endl;
		psi_momentum = FFT::FFT(psi);
		output_p << psi_momentum.transpose().format(CommaInitFmt) << std::endl;
	}
	std::cout << "Data generation finished. Input anything and press enter to exit...";

	output_x.close();
	output_p.close();
	char c;
	std::cin >> c;
	return 0;
}