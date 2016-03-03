#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/Geometry> 

namespace Eigen
{
	template<typename _Scalar>
class QuaternionAlias
{
	public:
		Quaternion<_Scalar> q;
		QuaternionAlias( void ) { }
		QuaternionAlias( _Scalar w, _Scalar x, _Scalar y, _Scalar z ) { q = Quaternion<_Scalar>( w, -x, -y, -z ); }
		QuaternionAlias( const Quaternion<_Scalar> & q_ ) { q = q_.conjugate(); }
		Quaternion<_Scalar> toQuat() { return q.conjugate(); }
		Quaternion<_Scalar>& getRaw() { return q; }
		void normalize() { q.normalize(); }
};

typedef QuaternionAlias<float> QuaternionAf;
typedef QuaternionAlias<double> QuaternionAd;

}

namespace ekf
{

class EstimatorPredictor
{
public:
	EstimatorPredictor();
	void PropagateState(const Eigen::Vector3d &omega_m, const Eigen::Vector3d &a_m);
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
			const Eigen::Vector3d &g_);
	void GetState(Eigen::Vector3d &p_i_w_,
			Eigen::Vector3d &v_i_w_,
			Eigen::Quaterniond &q_i_w_,
			Eigen::Vector3d &omega_i,
			Eigen::Vector3d &a_i);
	Eigen::Matrix<double,31,1> GetStateVector(void);


// protected:
	// States:
	Eigen::Vector3d p_i_w;
	Eigen::Vector3d v_i_w;
	Eigen::QuaternionAd q_i_w;
	Eigen::Vector3d b_omega;
	Eigen::Vector3d b_a;
	double lambda;
	Eigen::Vector3d p_c_i;
	Eigen::QuaternionAd q_c_i;
	Eigen::Vector3d p_w_v;
	Eigen::QuaternionAd q_w_v;

	// Calibration parameters
	double sq_sigma_omega;
	double sq_sigma_a;
	double sq_sigma_b_omega;
	double sq_sigma_b_a;
	double Delta_t;
	Eigen::Vector3d g;

	// Saved for trapetziodal integration
	Eigen::Vector3d omega_m_old;
	Eigen::Vector3d a_m_old;
};

class EstimatorFull : public EstimatorPredictor
{
public:
	EstimatorFull();

	// void PropagateState(const Eigen::Vector3d &omega_m, const Eigen::Vector3d &a_m); USE SUPERCLASS METHOD
	void SetCovarianceDiagonal(const Eigen::Matrix<double,28,1> &P_);
	void PropagateCovariance(const Eigen::Vector3d &omega_m, const Eigen::Vector3d &a_m);
	void UpdateCamera(const Eigen::Vector3d &p_c_v, const Eigen::Quaterniond &q_c_v,
		const Eigen::Matrix<double,6,6> &R,
		bool absolute, bool isKeyframe=false, double dLambda=0 );
	void UpdateKeyframe(void);
	Eigen::Matrix<double,28,1> GetCovarianceDiagonal(void);

//protected:
	// Covariance
	Eigen::Matrix<double,28,28> P;

//private:
	// Helper functions
	void ApplyCorrection(const Eigen::Matrix<double,28,1> &x);

	// Saved for calculating keyframe increments from absolute path
	Eigen::Vector3d p_kf_v;
	Eigen::QuaternionAd q_kf_v;

};

}