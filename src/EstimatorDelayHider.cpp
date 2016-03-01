#include "EstimatorDelayHider.h"
#include <iostream>


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
		const Eigen::Vector3d &p_w_v_,
		const Eigen::Quaterniond &q_w_v_)
{
	boost::mutex::scoped_lock lock(estimatorFullMutex);
	estimatorFull.SetState( p_i_w_, v_i_w_,	q_i_w_,	b_omega_,
			b_a_, lambda_, p_c_i_, q_c_i_, p_w_v_, q_w_v_ );
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
			const Eigen::Vector3d &a_m, double timeStamp)
{
	// Construct buffer object
	imuData_t data;
	data.omega_m = omega_m;
	data.a_m = a_m;
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
	std::cout << "Estimator thread started" << std::endl;
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
		while(true)
		{
			// Get latest imu
			imuData_t imuData;
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

			// Lock the estimator object and propagate
			boost::mutex::scoped_lock lockEstimator(estimatorFullMutex);
			estimatorFull.PropagateState( imuData.omega_m, imuData.a_m );
			estimatorFull.PropagateCovariance( imuData.omega_m, imuData.a_m );
			lockEstimator.unlock();
		}
		if (skip)
			continue;

		// Do update from camera
		boost::mutex::scoped_lock lockEstimator(estimatorFullMutex);
		estimatorFull.UpdateCamera(cameraData.p_c_v, cameraData.q_c_v, cameraData.R,
			absolute, cameraData.isNewKeyframe, Delta_lambda);
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

	std::cout << "Estimator thread ending" << std::endl;
}

void
EstimatorDelayHider::PredictorThread(void)
{
	std::cout << "Predictor thread started" << std::endl;
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


	std::cout << "Predictor thread ending" << std::endl;
}


}