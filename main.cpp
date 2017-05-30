#include "main.h"

int main() {

	double h_bar=1, mass=1;
	//data input
	std::cout << std::endl << "-------- Properties of space and time --------" << std::endl << std::endl;
	double step_t, step_x;
	INPUT("Input step in space: ", step_x);
	INPUT("Input step in time: ", step_t);
	int N_space = (int)(SYSTEM_SIZE/step_x); //number of space points
	int N_time = (int)(TIME_SIZE/step_t); //number of time points
	OUTPUT("Constant h bar: ",h_bar);
	OUTPUT("Constant mass: ",mass);
	OUTPUT("System spatial length: ",SYSTEM_SIZE);
	OUTPUT("System time length: ",TIME_SIZE);
	OUTPUT("Total steps in time: ",N_time);
	OUTPUT("Total steps in space: ",N_space);
	OUTPUT("Ratio step_x/step_t: ",step_x/step_t);
	
	VectorC space = utils::linspace(step_x, SYSTEM_SIZE);
	VectorC psi(N_space);

	std::cout << std::endl << "-------- Initial parameters of Gaussian wavepacket --------" << std::endl << std::endl;
	double momentum, initial_pos, spread;
	INPUT("Initial position: ", initial_pos);
	INPUT("Initial spread: ", spread);
	INPUT("Initial momentum: ", momentum);
	
	for(int i=0; i<N_space; i++) {
		dcplx val = utils::gauss_state(space(i), momentum, initial_pos, spread);
		psi(i) = val;
	}
	
	OUTPUT("Velocity: ", h_bar*momentum/mass);
	OUTPUT("Maximal free particle time window: ", h_bar*momentum/mass);
	
	std::cout << "-------- Potential function --------" << std::endl << std::endl;
	OUTPUT("Not yet implemented... ", " Sorry");
	VectorR pot = VectorR::Zero(N_space);

	
	//system init
	QMSystem system(psi, pot, step_x, step_t, h_bar, mass);
	
	//std::cout << system.get_state().cwiseAbs2() << std::endl;

	//output init
	std::ofstream output_x, output_p;
	output_x.open("out-space.csv");
	output_p.open("out-momentum.csv");
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ";", "", "", "", "", ""); //copied from reference
	
	output_x << psi.transpose().format(CommaInitFmt) << std::endl;
	VectorC psi_momentum = FFT::FFT(psi);
	output_p << psi_momentum.transpose().format(CommaInitFmt) << std::endl;
	for(int t=0; t<N_time;t++) {
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