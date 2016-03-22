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
			Eigen::Quaterniond(1,0.1,0,0),		// q_i_w
			Eigen::Vector3d(0,0,0),				// b_omega
			Eigen::Vector3d(0,0,0),				// b_a
			log(1),								// lambda
//			Eigen::Vector3d(0.045, 0, -0.05),	// p_c_i
//			Eigen::Quaterniond(0,1,1,0),		// q_c_i
			Eigen::Vector3d(0.1, 0, 0),	// p_c_i
			Eigen::Quaterniond(0,1,1,0),		// q_c_i
			Eigen::Vector3d(0,0,0),				// p_ikf_w
			Eigen::Quaterniond(1,0,0,0));		// q_ikf_w
	filter.SetCalibration(0.02*0.02,			// sq_sigma_omega
			0.05*0.05,							// sq_sigma_a
			0.0001*0.0001,						// sq_sigma_b_omega
			0.0001*0.0001,						// sq_sigma_a_omega
			1/400.0,							// Delta_t
			Eigen::Vector3d(0,0,9.82));			// g
	Eigen::Matrix<double,28,1> P;
	P << 0, 0, 0,								// p_i_w
	0.1, 0.1, 0.1,									// v_i_w
	0.3, 0.3, 0,								// q_i_w
	0.01, 0.01, 0.01,								// b_omega
	0.1, 0.1, 0.1,									// b_a
	log(2),											// lambda
	0.01, 0.01, 0.01,							// p_c_i
	0.01, 0.01, 0.01,								// q_c_i
	0, 0, 0,									// p_ikf_w
	0, 0, 0;									// q_ikf_w
	P = 1*P.cwiseProduct(P);
	filter.SetCovarianceDiagonal(P);
	filter.UpdateKeyframe();

	static Eigen::Vector3d p_kf_v(0,0,0);
	static Eigen::QuaternionAd q_kf_v(1,0,0,0);

	unsigned numSlamMeas = 0;
	for( std::string line; getline( slamFile, line ); )
	{
		std::stringstream ss;
		ss.str( line );

		Eigen::Vector3d p_c_v;
		Eigen::Quaterniond q_c_v;
		double slamTimeNs;	
		int newKf;

		char tmp;
		ss
			>> p_c_v(0) >> tmp >> p_c_v(1) >> tmp >> p_c_v(2) >> tmp
	    	>> q_c_v.w() >> tmp >> q_c_v.x() >> tmp >> q_c_v.y() >> tmp >> q_c_v.z() >> tmp
	    	>> newKf >> tmp
	    	>> slamTimeNs;

	    numSlamMeas++;
	    q_c_v.normalize(); // q is scaled to indicate scaling of keyframe. This might be extracted to set covariance
	    //q_c_v.conjugate();

	    Eigen::QuaternionAd q_c_kf(q_c_v.conjugate() * q_kf_v.toQuat().conjugate());
	    Eigen::Vector3d p_c_kf = filter.q_kf_v.toQuat().toRotationMatrix() * (p_c_v - filter.p_kf_v);
		

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
		    if (numSlamMeas>1)
		    {
		    	Eigen::Matrix3d C_q_ikf_w = filter.q_ikf_w.toQuat().toRotationMatrix();
				Eigen::Matrix3d C_q_i_w = filter.q_i_w.toQuat().toRotationMatrix();
				Eigen::Matrix3d C_q_c_i = filter.q_c_i.toQuat().toRotationMatrix();

				
		    	std::cout << filter.p_i_w.transpose() << " " << filter.q_i_w.q.coeffs().transpose() << " "
		   			<< filter.p_ikf_w.transpose() << " " << filter.q_ikf_w.q.coeffs().transpose() << " "
		    		<< filter.lambda << " " << exp(filter.lambda) << " "
		    		<< imuTimeNs << " "
		    		<< slamTimeNs << " "
		    		<< filter.GetStateVector().transpose() << " " << filter.GetCovarianceDiagonal().transpose() << " "
		    		<< p_c_kf.transpose() << " "
		    		<< (C_q_c_i*( C_q_ikf_w*( filter.p_i_w + C_q_i_w.transpose()*filter.p_c_i - filter.p_ikf_w ) - filter.p_c_i )*exp(filter.lambda)).transpose() << " "
		    		<< (filter.q_c_i.toQuat()*filter.q_i_w.toQuat()*filter.q_ikf_w.toQuat()).conjugate().coeffs().transpose() << " "
		    		<< std::endl;

		    		// << " " << filter.P.diagonal().transpose() << std::endl;
				filter.PropagateState( omega_m, a_m );
				filter.PropagateCovariance( omega_m, a_m );
			}

			if ((imuTimeNs+0.0/60.0)>=slamTimeNs)
				break;
		}

		Eigen::Matrix<double,6,6> R = R.Zero();
		R.diagonal()[0] = 0.01*0.01;
		R.diagonal()[1] = 0.01*0.01;
		R.diagonal()[2] = 0.01*0.01;
		R.diagonal()[3] = 0.1*0.1;
		R.diagonal()[4] = 0.1*0.1;
		R.diagonal()[5] = 0.01*0.01;
		if ((numSlamMeas%1 == 0)&&(numSlamMeas<60*10000))
		{

			filter.UpdateCamera(p_c_v, q_c_v, R, true, newKf, 0*0.001*0.001*log(2)*log(2));
			//filter.UpdateCamera(p_c_v, q_c_v, R, true, numSlamMeas%50 == 0, 0*0.001*0.001*log(2)*log(2));
		}
		//filter.UpdateCamera(Eigen::Vector3d(0,0,0), Eigen::Quaterniond(1,0,0,0), R);

	}
	return 0;
}