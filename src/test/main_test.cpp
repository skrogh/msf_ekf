#include "EstimatorBase.h"
#include "EstimatorDelayHider.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/thread.hpp>

int
 main(int argc, char *argv[])
{
	std::ifstream testFile;
	if (argc > 1)
	{
		testFile.open( argv[1] );
	}
	else
	{
		std::cerr << "Please provide file for testing!" << std::endl;
		exit(1);
	}
	if (!testFile.is_open())
	{
		std::cerr << "Error opening file: " << argv[1] << std::endl;
		exit(1);
	}

	ekf::EstimatorFull filter;

	unsigned i = 0;
	for( std::string line; getline( testFile, line ); )
	{
		std::stringstream ss;
		ss.str( line );

		char tmp;
		
		Eigen::Vector3d a_m;
		Eigen::Vector3d omega_m;
		Eigen::Vector3d p_c_v;
		Eigen::Quaterniond q_c_v;

	    ss >> omega_m(0) >> tmp >> omega_m(1) >> tmp >> omega_m(2) >> tmp
	    >> a_m(0) >> tmp >> a_m(1) >> tmp >> a_m(2) >> tmp 
	    >> p_c_v(0) >> tmp >> p_c_v(1) >> tmp >> p_c_v(2) >> tmp
	    >> q_c_v.coeffs()(3) >> tmp >> q_c_v.coeffs()(0) >> tmp >> q_c_v.coeffs()(1) >> tmp >> q_c_v.coeffs()(2);

		//std::cout << "omega_m: \n" << omega_m << std::endl;
		//std::cout << "a_m: \n" << a_m << std::endl;
		std::cout << filter.p_i_w.transpose() << " " << filter.q_i_w.coeffs().transpose() << " " << filter.lambda << " " << exp(filter.lambda) << " " << filter.P.diagonal().transpose() << std::endl;


		filter.PropagateState( omega_m, a_m );
		filter.PropagateCovariance( omega_m, a_m );

		if ( i%(400/60) == 1)
		{
			Eigen::Matrix<double,6,6> R = R.Zero();
			R.diagonal()[0] = 0.01;
			R.diagonal()[1] = 0.01;
			R.diagonal()[2] = 0.01;
			R.diagonal()[3] = 0.1;
			R.diagonal()[4] = 0.1;
			R.diagonal()[5] = 0.1;

			filter.UpdateCamera(p_c_v, q_c_v, R);
		}

		i++;
	}
	return 0;
}