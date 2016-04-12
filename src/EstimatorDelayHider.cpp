#include "EstimatorDelayHider.h"
#include <iostream>
#include <boost/circular_buffer.hpp>

#include <unistd.h>
#include <sys/syscall.h>

// Helper function for signum
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double calcLambda(double S_xx, double S_yy, double S_xy,
	double sigma_x, double sigma_y)
{
	S_xx *= sigma_y*sigma_y;
	S_yy *= sigma_x*sigma_x;
	S_xy *= sigma_x*sigma_y;
	return (S_xx - S_yy + sgn(S_xy)*sqrt( (S_xx - S_yy)*(S_xx - S_yy) + 4*S_xy*S_xy)) / (2/sigma_x*sigma_y*S_xy);
}

namespace ekf
{

EstimatorDelayHider::EstimatorDelayHider()
{

}

EstimatorDelayHider::~EstimatorDelayHider()
{
	boost::mutex::scoped_lock lock(threadsShouldExitMutex);
	threadsShouldExit = true;
	lock.unlock();
	cameraData_t cameraDummy;
	imuData_t imuDummy;
	imuPredictorBuffer.push(imuDummy);
	cameraInBuffer.push(cameraDummy);
	estimatorThread->join();
	predictorThread->join();
	delete estimatorThread;
	delete predictorThread;
}

void
EstimatorDelayHider::SetCovarianceDiagonal(const Eigen::Matrix<double,28,1> &P_)
{
	boost::mutex::scoped_lock lock(estimatorFullMutex);
	estimatorFull.SetCovarianceDiagonal( P_ );
}

void
EstimatorDelayHider::SetState(const Eigen::Vector3d &p_i_w_,
		const Eigen::Vector3d &v_i_w_,
		const Eigen::Quaterniond &q_i_w_,
		const Eigen::Vector3d &b_omega_,
		const Eigen::Vector3d &b_a_,
		double lambda_,
		const Eigen::Vector3d &p_c_i_,
		const Eigen::Quaterniond &q_c_i_,
		const Eigen::Vector3d &p_ikf_w_,
		const Eigen::Quaterniond &q_ikf_w_)
{
	boost::mutex::scoped_lock lock(estimatorFullMutex);
	estimatorFull.SetState( p_i_w_, v_i_w_,	q_i_w_,	b_omega_,
			b_a_, lambda_, p_c_i_, q_c_i_, p_ikf_w_, q_ikf_w_ );
}

void
EstimatorDelayHider::SetCalibration(double sq_sigma_omega_,
		double sq_sigma_a_,
		double sq_sigma_b_omega_,
		double sq_sigma_b_a_,
		double Delta_t_,
		const Eigen::Vector3d &g_,
		double Delta_lambda_,
		bool absolute_)
{
	Delta_lambda = Delta_lambda_;
	absolute = absolute_;
	boost::mutex::scoped_lock lock(estimatorFullMutex);
	estimatorFull.SetCalibration(sq_sigma_omega_, sq_sigma_a_,
			sq_sigma_b_omega_, sq_sigma_b_a_, Delta_t_, g_);
}

void
EstimatorDelayHider::CameraMeasurement(const Eigen::Vector3d &p_c_v,
		const Eigen::Quaterniond &q_c_v,
		const Eigen::Matrix<double,6,6> &R,
		bool isNewKeyframe,
		double timeStamp)
{
	cameraData_t data;
	data.p_c_v = p_c_v;	
	data.q_c_v = q_c_v;
	data.R = R;
	data.isNewKeyframe = isNewKeyframe;
	data.timeStamp = timeStamp;

	cameraInBuffer.push(data);
}

void
EstimatorDelayHider::ImuMeasurement(const Eigen::Vector3d &omega_m,
			const Eigen::Vector3d &a_m, double dist, bool distValid, double timeStamp)
{
	// Construct buffer object
	imuData_t data;
	data.omega_m = omega_m;
	data.a_m = a_m;
	data.dist = dist;
	data.distValid = distValid;
	data.timeStamp = timeStamp;

	// Check if predictor finished
	// Check is done by checking if the predictor is running and has no more steps to process
	// Replace current estimator if it did
	boost::mutex::scoped_lock lock(estimatorPredictorMutex);
	if ( predictorRunning && imuPredictorBuffer.empty() )
	{
		estimatorPredictorCurrent = estimatorPredictor;
		predictorRunning = false;
	}
	else if ( predictorRunning )
	{
		imuPredictorBuffer.push(data);
	}
	// Push data to updater thread. Do it inside predictor thread lock, as this prevents a new copy to the predictor, while we add elements to the buffers.
	imuInBuffer.push(data);
	lock.unlock();

	// Propagate update for main thread (the one that is actually seen and used)
	estimatorPredictorCurrent.PropagateState( omega_m, a_m );
}

void 
EstimatorDelayHider::Start(void)
{
	estimatorThread = new boost::thread(boost::bind(&EstimatorDelayHider::EstimatorThread, this)); // Spawn thread running the estimatorThread
	predictorThread = new boost::thread(boost::bind(&EstimatorDelayHider::PredictorThread, this)); // Spawn thread running the predictorThread
}

void
EstimatorDelayHider::GetState(Eigen::Vector3d &p_i_w_,
			Eigen::Vector3d &v_i_w_,
			Eigen::Quaterniond &q_i_w_,
			Eigen::Vector3d &omega_i,
			Eigen::Vector3d &a_i)
{
	estimatorPredictorCurrent.GetState(p_i_w_, v_i_w_, q_i_w_, omega_i, a_i);
}

void
EstimatorDelayHider::EstimatorThread(void)
{
	bool initialized = false;
	std::clog << "Estimator thread started. PID: " << (int) syscall(SYS_gettid) << std::endl;

	// Local variables for scale estimation
	int numValidDist = 0;
	double distanceBuffer[50];

	int windowSize = 10;
	boost::circular_buffer<double> distanceRing(windowSize);
	boost::circular_buffer<Eigen::Vector3d> p_c_vRing(windowSize);
	boost::circular_buffer<Eigen::QuaternionAd> q_w_vRing(windowSize);
	double S_xx = 0;
	double S_yy = 0;
	double S_xy = 0;

	double sigma_x = 1;
	double sigma_y = 1;
	double tau = 1 - 1/200;	

	while ( 1 )
	{
		boost::mutex::scoped_lock lockExit(threadsShouldExitMutex);
		if( threadsShouldExit )
			break;
		lockExit.unlock();

		cameraData_t cameraData;
		cameraInBuffer.popBlocking(cameraData);

		// Propagate up to new measurement
		bool skip = false;
		imuData_t imuData;
		while(true)
		{
			// Get latest imu
			if (!imuInBuffer.peek(imuData))
			{
				std::cerr << "ESTIMATOR THREAD ERROR! IMU buffer empty after getting new camera image, this should not happen!" << std::endl
				<< "Will skip camera image for update" << std::endl;
				skip = true;
				break;
			}
			// Check if imu data is newer than camera. If it is, we're done propagating
			if (imuData.timeStamp>cameraData.timeStamp) // TODO: check if times are valid (if they are not skip the update)
				break;

			// Imu measurement was older than camera. Grab it and continue with prediction
			if (!imuInBuffer.popNonBlocking(imuData))
			{
				std::cerr << "ESTIMATOR THREAD ERROR! IMU buffer empty after check for not empty!" << std::endl
				<< "Someone else is consuming!" << std::endl;
				skip = true;
				break;
			}

			if (initialized) //skip untill first measurment
			{
				// Lock the estimator object and propagate
				boost::mutex::scoped_lock lockEstimator(estimatorFullMutex);
				if ( imuData.distValid )
				{
					estimatorFull.UpdateAltitude( imuData.dist, 0.01*0.01 );
				}
				estimatorFull.PropagateState( imuData.omega_m, imuData.a_m );
				estimatorFull.PropagateCovariance( imuData.omega_m, imuData.a_m );
				lockEstimator.unlock();
				if (imuData.distValid)
					distanceBuffer[numValidDist++] = imuData.dist;
			}
		}
		if (skip)
			continue;

		initialized = true;

		// Do update from camera
		boost::mutex::scoped_lock lockEstimator(estimatorFullMutex);
		estimatorFull.UpdateCamera(cameraData.p_c_v, cameraData.q_c_v, cameraData.R,
			absolute, cameraData.isNewKeyframe, Delta_lambda);
		
		static int i = 0;
		if (i++%1==0)
		{
			Eigen::QuaternionAd q_c_kf(cameraData.q_c_v.conjugate() * estimatorFull.q_kf_v.toQuat().conjugate());
    		Eigen::Vector3d p_c_kf = estimatorFull.q_kf_v.toQuat().toRotationMatrix() * (cameraData.p_c_v - estimatorFull.p_kf_v);
	    	Eigen::Matrix3d C_q_w_v = estimatorFull.q_ikf_w.toQuat().toRotationMatrix();
			Eigen::Matrix3d C_q_i_w = estimatorFull.q_i_w.toQuat().toRotationMatrix();
			Eigen::Matrix3d C_q_c_i = estimatorFull.q_c_i.toQuat().toRotationMatrix();

	    	std::cout << estimatorFull.p_i_w.transpose() << " " << estimatorFull.q_i_w.q.coeffs().transpose() << " "
	   			<< estimatorFull.p_ikf_w.transpose() << " " << estimatorFull.q_ikf_w.q.coeffs().transpose() << " "
	    		<< estimatorFull.lambda << " " << exp(estimatorFull.lambda) << " "
	    		<< imuData.timeStamp << " "
	    		<< cameraData.timeStamp << " "
	    		<< estimatorFull.GetStateVector().transpose() << " " << estimatorFull.GetCovarianceDiagonal().transpose() << " "
	    		<< p_c_kf.transpose() << " " << ( C_q_w_v.transpose()*(estimatorFull.p_i_w + C_q_i_w.transpose()*estimatorFull.p_c_i) - estimatorFull.p_ikf_w ).transpose() * exp(estimatorFull.lambda) << " "
	    		<< (estimatorFull.q_c_i.toQuat()*estimatorFull.q_i_w.toQuat()*estimatorFull.q_ikf_w.toQuat()).conjugate().coeffs().transpose() << " "
	    		<< cameraData.p_c_v.transpose() << " "
	    		<< std::endl;

		}

	   	Eigen::QuaternionAd q_w_v( ( estimatorFull.q_c_i.toQuat() * estimatorFull.q_i_w.toQuat() ).conjugate() * cameraData.q_c_v.conjugate() );
		Eigen::Vector3d p_c_w = q_w_v.toQuat().toRotationMatrix() * cameraData.p_c_v;

		if (numValidDist)
		{
			distanceRing.push_back(distanceBuffer[numValidDist-1]);
			p_c_vRing.push_back(cameraData.p_c_v);
			q_w_vRing.push_back(q_w_v);
		}

	    if (distanceRing.full()&& numValidDist)
	    {
	    	double x_ = distanceRing.back() - distanceRing.front();
	    	Eigen::Vector3d y_ =  q_w_vRing.back().toQuat().toRotationMatrix()*(p_c_vRing.back() - p_c_vRing.front());

	    	S_xx = S_xx*tau + x_*x_;
	    	S_yy = S_yy*tau + y_(2)*y_(2);
	    	S_xy = S_xy*tau + x_*y_(2);
	    	

	    	double lambda = calcLambda(S_xx, S_yy, S_xy, sigma_x, sigma_y);

	    	// Apply lambda "measurement"
	   		// estimatorFull.lambda = lambda;

	    }
	    numValidDist = 0;


		lockEstimator.unlock();


		// Push newly updated state into predictor to catch up to current non-delayed estimate
		boost::mutex::scoped_lock lockPredictor(estimatorPredictorMutex);
		if( predictorRunning )
		{
			std::cerr << "ESTIMATOR THREAD ERROR! Predictor is too slow. Did not finish before new measurement!" << std::endl;
			continue;
		}
		estimatorPredictor = estimatorFull;
		imuPredictorBuffer = imuInBuffer;
		predictorRunning = true;
		lockPredictor.unlock();
	}

	std::clog << "Estimator thread ending" << std::endl;
}

void
EstimatorDelayHider::PredictorThread(void)
{
	std::clog << "Predictor thread started. PID: " << (int) syscall(SYS_gettid) << std::endl;
	while ( 1 )
	{
		boost::mutex::scoped_lock lockExit(threadsShouldExitMutex);
		if( threadsShouldExit )
			break;
		lockExit.unlock();

		imuData_t data;
		imuPredictorBuffer.blockWhileEmpty();

		boost::mutex::scoped_lock lock(estimatorPredictorMutex);
		if( !imuPredictorBuffer.popNonBlocking(data) )
			std::cerr << "PREDICTOR THREAD ERROR! Someone else is consuming the predictor catch-up buffer!" << std::endl;
		estimatorPredictor.PropagateState( data.omega_m, data.a_m );
		lock.unlock();

	}


	std::clog << "Predictor thread ending" << std::endl;
}


}