#include "EstimatorBase.h"
#include "EstimatorDelayHider.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/thread.hpp>

int
 main(int argc, char *argv[])
{
	std::ifstream imuFile;
	std::ifstream slamFile;
	if (argc > 2)
	{
		imuFile.open( argv[1] );
		slamFile.open( argv[2] );
	}
	else
	{
		std::cerr << "Please provide files for testing!" << std::endl
		<< "Format is: ekf_test <imuFile> <slamFile>" << std::endl;
		exit(1);
	}
	if (!imuFile.is_open())
	{
		std::cerr << "Error opening file: " << argv[1] << std::endl;
		exit(1);
	}
	if (!slamFile.is_open())
	{
		std::cerr << "Error opening file: " << argv[2] << std::endl;
		exit(1);
	}

	ekf::EstimatorFull filter;
	filter.SetState( Eigen::Vector3d(0, 0, 0), 	// p_i_w
			Eigen::Vector3d(0,0,0),				// v_i_w
			Eigen::Quaterniond(1,0,0,0),		// q_i_w
			Eigen::Vector3d(0,0,0),				// b_omega
			Eigen::Vector3d(0,0,0),				// b_a
			log(1),								// lambda
			Eigen::Vector3d(0,0.1,0),// p_c_i
			Eigen::Quaterniond(1,0,0,0),		// q_c_i
			Eigen::Vector3d(0,-0.1,0),				// p_w_v
			Eigen::Quaterniond(1,0,0,0));		// q_w_v
	filter.SetCalibration(0.0002*400.0,			// sq_sigma_omega
			0.005*400.0,							// sq_sigma_a
			0.002/400,							// sq_sigma_b_omega
			0.05/400,							// sq_sigma_a_omega
			1/400.0,							// Delta_t
			Eigen::Vector3d(0,0,9.82));			// g
	Eigen::Matrix<double,28,1> P;
	P << 0, 0, 0,								// p_i_w
	0, 0, 0,									// v_i_w
	0.5, 0.5, 0,								// q_i_w
	0, 0, 0,								// b_omega
	0, 0, 0,									// b_a
	0,											// lambda
	0, 0, 0,							// p_c_i
	0, 0, 0,								// q_c_i
	0, 0, 0,									// p_w_v
	0, 0, 0;									// q_w_v
	P = 9.0*P.cwiseProduct(P);
	filter.SetCovarianceDiagonal(P);
	//filter.UpdateKeyframe();

	unsigned numSlamMeas = 0;
	for( std::string line; getline( slamFile, line ); )
	{
		std::stringstream ss;
		ss.str( line );

		Eigen::Vector3d p_c_v;
		Eigen::Quaterniond q_c_v;
		double slamTimeNs;		

		char tmp;
		ss
			>> p_c_v(0) >> tmp >> p_c_v(1) >> tmp >> p_c_v(2) >> tmp
	    	>> q_c_v.coeffs()(3) >> tmp >> q_c_v.coeffs()(0) >> tmp >> q_c_v.coeffs()(1) >> tmp >> q_c_v.coeffs()(2) >> tmp
	    	>> slamTimeNs;
	    numSlamMeas++;
	    q_c_v.normalize(); // q is scaled to indicate scaling of keyframe. This might be extracted to set covariance

		// propagate
		for( std::string imuLine; getline( imuFile, imuLine ); )
		{
			std::stringstream ss;
			ss.str( imuLine );

			Eigen::Vector3d a_m;
			Eigen::Vector3d omega_m;
			double imuTimeNs;		

			char tmp;
			ss
				>> omega_m(0) >> tmp >> omega_m(1) >> tmp >> omega_m(2) >> tmp
	    		>> a_m(0) >> tmp >> a_m(1) >> tmp >> a_m(2) >> tmp 
		    	>> imuTimeNs;
		    if (numSlamMeas>0)
		    {
		    	Eigen::Matrix3d C_q_w_v = filter.q_w_v.toRotationMatrix();
				Eigen::Matrix3d C_q_i_w = filter.q_i_w.toRotationMatrix();
				Eigen::Matrix3d C_q_c_i = filter.q_c_i.toRotationMatrix();

		    	std::cout << filter.p_i_w.transpose() << " " << filter.q_i_w.coeffs().transpose() << " "
		   			<< filter.p_w_v.transpose() << " " << filter.q_w_v.coeffs().transpose() << " "
		    		<< filter.lambda << " " << exp(filter.lambda) << " "
		    		<< imuTimeNs << " "
		    		<< slamTimeNs << " "
		    		<< filter.GetStateVector().transpose() << " " << filter.GetCovarianceDiagonal().transpose() << " "
		    		<< p_c_v.transpose() << " " << ( C_q_w_v.transpose()*(filter.p_i_w + C_q_i_w.transpose()*filter.p_c_i) + filter.p_w_v ).transpose() * exp(filter.lambda) << " "
		    		<< std::endl;

		    		// << " " << filter.P.diagonal().transpose() << std::endl;
				filter.PropagateState( omega_m, a_m );
				filter.PropagateCovariance( omega_m, a_m );
			}

			if ((imuTimeNs)>=slamTimeNs)
				break;
		}

		Eigen::Matrix<double,6,6> R = R.Zero();
		R.diagonal()[0] = 0.001;
		R.diagonal()[1] = 0.001;
		R.diagonal()[2] = 0.001;
		R.diagonal()[3] = 1000.1;
		R.diagonal()[4] = 1000.1;
		R.diagonal()[5] = 1000.1;
		if (numSlamMeas%2 == 0)
		{
			//filter.UpdateCamera(p_c_v, q_c_v, R);
			//if (numSlamMeas/10>20)
				//filter.UpdateKeyframe();
			filter.P(15,15) += 0; // lambda
			filter.P(22,22) += 0;
			filter.P(23,23) += 0;
			filter.P(24,24) += 0;
			
			filter.P(25,25) += 0;
			filter.P(26,26) += 0;
			filter.P(27,27) += 0;
		}
		//filter.UpdateCamera(Eigen::Vector3d(0,0,0), Eigen::Quaterniond(1,0,0,0), R);

	}
	return 0;
}