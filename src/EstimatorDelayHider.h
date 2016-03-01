#pragma once
#include "EstimatorBase.h"
#include "ConcurrentBuffers.h"
#include <boost/thread.hpp>



namespace ekf
{

struct imuData_t
{
	Eigen::Vector3d omega_m;	// Measured angular velocity in rad/s
	Eigen::Vector3d a_m;		// Measured acceleration in m/s²
	double timeStamp; 				// Time of measurement
};

struct cameraData_t
{
	Eigen::Vector3d p_c_v;			// Measured position in visual frame
	Eigen::Quaterniond q_c_v;		// Measured Quaternion in visual frame
	Eigen::Matrix<double,6,6> R;	// Measurement covariance
	bool isNewKeyframe;				// Set if measurement is a new keyframe
	double timeStamp; 	// Time of measurement
};


class EstimatorDelayHider
{
public:
	EstimatorDelayHider();
	~EstimatorDelayHider();

	// Warning! Setting state/covariance/calibration while filter is running might have undesired effects!
	void SetCovarianceDiagonal(const Eigen::Matrix<double,28,1> &P_);
	void SetState(const Eigen::Vector3d &p_i_w_,
			const Eigen::Vector3d &v_i_w_,
			const Eigen::Quaterniond &q_i_w_,
			const Eigen::Vector3d &b_omega_,
			const Eigen::Vector3d &b_a_,
			double lambda_,
			const Eigen::Vector3d &p_c_i_,
			const Eigen::Quaterniond &q_c_i_,
			const Eigen::Vector3d &p_w_v_,
			const Eigen::Quaterniond &q_w_v_);
	void SetCalibration(double sq_sigma_omega_,
			double sq_sigma_a_,
			double sq_sigma_b_omega_,
			double sq_sigma_b_a_,
			double Delta_t_,
			const Eigen::Vector3d &g_,
			double Delta_lambda_,
			bool absolute_);
	void Start(void);
	void GetState(Eigen::Vector3d &p_i_w_,
			Eigen::Vector3d &v_i_w_,
			Eigen::Quaterniond &q_i_w_,
			Eigen::Vector3d &omega_i,
			Eigen::Vector3d &a_i);
	
	void CameraMeasurement(const Eigen::Vector3d &p_c_v,
			const Eigen::Quaterniond &q_c_v,
			const Eigen::Matrix<double,6,6> &R,
			bool isNewKeyframe,
			double timeStamp);

	void ImuMeasurement(const Eigen::Vector3d &omega_m,
			const Eigen::Vector3d &a_m, double timeStamp); // Return state somehow?


//protected:
	// Calibration values
	double Delta_lambda;
	bool absolute;

	// Buffers
	ConcurrentQueue<imuData_t> imuInBuffer;
	ConcurrentQueue<cameraData_t> cameraInBuffer;
	ConcurrentQueue<imuData_t> imuPredictorBuffer;

	// Estimators
	EstimatorPredictor estimatorPredictor;			// Used for "catching up" the predicted state, when new iamges arrive
	EstimatorPredictor estimatorPredictorCurrent;	// New measurements are estimated in this imediately (not in a thread, but called from "ImuMeasurement function")
	EstimatorFull estimatorFull;					// Full state predictor with covariance. Run in a thread.

	// Threads
	void EstimatorThread(void);
	void PredictorThread(void);
	boost::thread* estimatorThread;
	boost::thread* predictorThread;
	// Threads internal variables
	bool predictorRunning = false;


	// Mutexes
	mutable boost::mutex estimatorFullMutex;
	mutable boost::mutex estimatorPredictorMutex;
	mutable boost::mutex threadsShouldExitMutex;

	// Exit thread variable
	bool threadsShouldExit = false; //todo wrap in mutex?

};




}