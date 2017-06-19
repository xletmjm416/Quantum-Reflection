#include "main.h"

int main() {
	
	//units: eV, q_e, fs, microm
	double h_bar=0.6582, mass=5.685;
	//data input
	std::cout << std::endl << "-------- Properties of space and time --------" << std::endl << std::endl;
	double step_t, step_x, size_x, size_t;
	int N_space, N_time;
	INPUT("Input total size in space (nm): ", size_x);
	INPUT("Input total size in time (fs): ", size_t);
	INPUT("Input total steps in space (nm): ", N_space);
	INPUT("Input total steps in time (fs): ", N_time);
	step_x = (double)size_x/(double)N_space;
	step_t = (double)size_t/(double)N_time;
	OUTPUT("Constant h bar (units eV / fs): ",h_bar);
	OUTPUT("Constant mass (units m_e): ",mass);
	OUTPUT("System spatial length (nm): ",size_x);
	OUTPUT("System time length (fs): ",size_t);
	OUTPUT("Total steps in time: ",N_time);
	OUTPUT("Total steps in space: ",N_space);
	OUTPUT("Step in space (nm): ", step_x);
	OUTPUT("Step in time (fs): ", step_t);
	OUTPUT("Ratio step_x/step_t (nm/fs): ",step_x/step_t);
	OUTPUT("Critical wavenumber (wavelength = system size) (1/nm): ",2*PI_CONST/size_x);
	OUTPUT("Nyquist wavenumber (wavelength = space step) (1/nm): ",PI_CONST/step_x);
	
	VectorC space = utils::linspace(step_x, size_x);
	VectorC psi(N_space);

	std::cout << std::endl << "-------- Initial parameters of Gaussian wavepacket --------" << std::endl << std::endl;
	double initial_pos, wavenumber, spread;
	INPUT("Initial position (nm): ", initial_pos);
	INPUT("Initial spread (nm): ", spread);
	INPUT("Initial wavenumber (1/nm): ", wavenumber);
	
	for(int i=0; i<N_space; i++) {
		dcplx val = utils::gauss_function(space(i), wavenumber, initial_pos, spread);
		//dcplx val = utils::step_function(space(i), 50, 20);
		psi(i) = val;
	}
	
	OUTPUT("Momentum (mass nm / fs): ", h_bar*wavenumber);
	OUTPUT("Momentum spread (mass nm / fs): ", h_bar/(2*spread));
	OUTPUT("Velocity (nm / fs): ", h_bar*wavenumber/mass);
	OUTPUT("Velocity spread (mass nm / fs): ", h_bar/(2*mass*spread));
	double energy_expectation = h_bar*h_bar/(2*mass)*wavenumber*wavenumber;
	OUTPUT("Energy ex. value (eV): ", energy_expectation);
	OUTPUT("Energy spread (eV): ", h_bar/(2*spread)*h_bar/(2*spread)/(2*mass));
	
	std::cout << std::endl << "-------- Potential function --------" << std::endl << std::endl;
	double height; double centre; short orient;
	INPUT("Input height of the barrier (eV): ", height);
	INPUT("Input centre of the barrier (nm): ", centre);
	INPUT("Input orientation (0 - left higher, 1- left lower): ", orient);
	VectorR pot(N_space);
	for(int i=0; i<N_space; i++) {
		dcplx val = utils::heaviside_function(space(i), centre, height, orient);
		pot(i) = val.real();
	}
	if (energy_expectation > height) {
	double other_wavenumber = std::pow(2*mass*(energy_expectation-height),0.5)/h_bar;
	OUTPUT("Wavenumber on top of barrier: ", other_wavenumber);
	OUTPUT("Reflection: ", 4*wavenumber*other_wavenumber/((wavenumber+other_wavenumber)*(wavenumber+other_wavenumber)));
	OUTPUT("Transmission: ", 1- 4*wavenumber*other_wavenumber/((wavenumber+other_wavenumber)*(wavenumber+other_wavenumber)));
	}
	else {
	OUTPUT("Decay constant: ", std::pow(2*mass*(height-energy_expectation),0.5)/h_bar);
	}
	
	//system init
	QMSystem system(psi, pot, step_x, step_t, h_bar, mass, size_x, size_t);
	
	//output init
	std::ofstream output_x, output_p;
	output_x.open("out-space.csv");
	//output_p.open("out-momentum.csv");
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ";", "", "", "", "", ""); //copied from reference
	
	output_x << h_bar << "," << mass << "," << step_x << "," << step_t << "," << size_x << "," << size_t << "," << N_space << "," << N_time << "," << initial_pos << "," << wavenumber << "," << spread << std::endl;
	//output_p << h_bar << "," << mass << "," << step_x << "," << step_t << "," << size_x << "," << size_t << "," << initial_pos << "," << wavenumber << "," << spread << std::endl;
	output_x << psi.transpose().format(CommaInitFmt) << std::endl;
	//VectorC psi_momentum = FFT::DITFFT(psi);
	//output_p << psi_momentum.transpose().format(CommaInitFmt) << std::endl;
	for(int t=0; t<N_time-1;t++) {
		psi = system.cranknicolson();
		output_x << psi.transpose().format(CommaInitFmt) << std::endl;
		//psi_momentum = FFT::DITFFT(psi);
		//output_p << psi_momentum.transpose().format(CommaInitFmt) << std::endl;
	}
	std::cout << "Data generation finished. Input anything and press enter to exit... ";

	output_x.close();
	//output_p.close();
	
	char c;
	std::cin >> c;
	return 0;
}