#include "EstimatorBase.h"
#include <iostream>
#include <math.h>       /* exp */


namespace Eigen
{
	Eigen::Matrix3d
	crossMat(const Eigen::Vector3d &t)
	{
		Eigen::Matrix3d t_hat;
		t_hat << 0, -t(2), t(1),
			t(2), 0, -t(0),
			-t(1), t(0), 0;
		return t_hat; 
	}
	Eigen::Matrix4d
	Omega( const Eigen::Vector3d& v ){
	Eigen::Matrix4d m;
	m <<     0,  v(2), -v(1),  v(0),
	     -v(2),     0,  v(0),  v(1),
	      v(1), -v(0),     0,  v(2),
	     -v(0), -v(1), -v(2),     0;
	return m;
}
}

using namespace Eigen;

namespace ekf
{

EstimatorPredictor::EstimatorPredictor()
{
	// States:
	p_i_w = Vector3d(0, 0, 0);
	v_i_w = Vector3d(0, 0, 0);
	q_i_w = QuaternionAd(1, 0, 0, 0);
	b_omega = Vector3d(0, 0, 0);
	b_a = Vector3d(0, 0, 0);
	lambda = log(1); // TOD: evaluate if lin or logspace is best
	p_c_i = Vector3d(0, 0, 0);
	q_c_i = QuaternionAd(1, 0, 0, 0);
	p_ikf_w = Vector3d(0, 0, 0);
	q_ikf_w = QuaternionAd(1, 0, 0, 0);

	// Calibration parameters
	g << 0, 0, 9.82;
	sq_sigma_omega = 1;
	sq_sigma_a = 1;
	sq_sigma_b_omega = 0.1;
	sq_sigma_b_a = 0.1;
	Delta_t = 1/400.0;

	// Saved for trapetziodal integration
	omega_m_old << 0, 0, 0;
	a_m_old = g;
	v_i_w_old = v_i_w;
	p_i_w_old = p_i_w;
	q_i_w_old = q_i_w;
}

void
EstimatorPredictor::SetState(const Eigen::Vector3d &p_i_w_,
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
	p_i_w = p_i_w_;
	v_i_w = v_i_w_;
	q_i_w = Eigen::QuaternionAd( q_i_w_ );
	q_i_w.normalize();
	b_omega = b_omega_;
	b_a = b_a_;
	lambda = lambda_;
	p_c_i = p_c_i_;
	q_c_i = Eigen::QuaternionAd( q_c_i_ );
	q_c_i.normalize();
	p_ikf_w = p_ikf_w_;
	q_ikf_w = Eigen::QuaternionAd( q_ikf_w_ );
	q_ikf_w.normalize();
}

void
EstimatorPredictor::SetCalibration(double sq_sigma_omega_,
		double sq_sigma_a_,
		double sq_sigma_b_omega_,
		double sq_sigma_b_a_,
		double Delta_t_,
		const Eigen::Vector3d &g_)
{
	sq_sigma_omega = sq_sigma_omega_;
	sq_sigma_a = sq_sigma_a_;
	sq_sigma_b_omega = sq_sigma_b_omega_;
	sq_sigma_b_a = sq_sigma_b_a_;
	Delta_t = Delta_t_;
	g = g_;
}

void
EstimatorPredictor::PropagateState(const Eigen::Vector3d &omega_m, const Eigen::Vector3d &a_m)
{
	Matrix3d C_q_i_w_old = q_i_w.toQuat().toRotationMatrix();
	v_i_w_old = v_i_w;
	p_i_w_old = p_i_w;
	q_i_w_old = q_i_w;

	Vector3d a_hat_old = a_m_old - b_a;
	Vector3d a_hat = a_m - b_a;
	Vector3d omega_hat_old = omega_m_old - b_omega;
	Vector3d omega_hat = omega_m - b_omega;
/*
	// (eye(4) - Delta_t/4*Omega(omega_m - b_omega))*q_i_w_new == (eye(4) + Delta_t/4*Omega( omega_m_old - b_omega))*q_i_w_old
	Vector3d tmp = Delta_t*(omega_m_old - b_omega)/2.0;
	Quaterniond Omega( 0, tmp(0), tmp(1), tmp(2) );
	q_i_w.q.coeffs() += (Omega.conjugate()*q_i_w.toQuat()).conjugate().coeffs();
	q_i_w.normalize();

	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();

	v_i_w += ( C_q_i_w_old.transpose()*(a_m_old - b_a) - g + C_q_i_w.transpose()*(a_m - b_a) - g )*Delta_t/2.0;
	p_i_w += (v_i_w_old + v_i_w)*Delta_t/2.0;
*/
	// Rotation
	Vector4d q0 = Vector4d( 0, 0, 0, 1 );
	Vector4d k1 = Omega( omega_hat_old ) * q0 / 2.0;
	Vector4d k2 = Omega( ( omega_hat_old + omega_hat ) / 2.0 ) * ( q0 + Delta_t/2.0 * k1 ) / 2.0;
	Vector4d k3 = Omega( ( omega_hat_old + omega_hat ) / 2.0 ) * ( q0 + Delta_t/2.0 * k2 ) / 2.0;
	Vector4d k4 = Omega( omega_hat ) * ( q0 + Delta_t * k3 ) / 2.0;

	Quaterniond q_i1_i(
			Vector4d( 0, 0, 0, 1 )
			+ Delta_t/6.0 * ( k1 + 2*k2 + 2*k3 + k4 )
	);
	Matrix3d C_i1_i = q_i1_i.conjugate().toRotationMatrix();

	q_i_w = QuaternionAd( q_i1_i.conjugate() * q_i_w.toQuat() );
	q_i_w.normalize();
	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();

	// Translation
	Vector3d a_hat_w = C_q_i_w.transpose()*a_hat - g;

	Vector3d s = Delta_t/2.0 * (
			C_i1_i.transpose()*a_hat + a_hat_old
	);

	Vector3d y = Delta_t/2.0 * s;

	v_i_w = v_i_w_old + C_q_i_w.transpose()*s - g * Delta_t;

	p_i_w = p_i_w_old + v_i_w_old * Delta_t
			+ C_q_i_w.transpose()*y
			- g * Delta_t*Delta_t/2.0;

	omega_m_old = omega_m;
	a_m_old = a_m;
}

void
EstimatorPredictor::GetState(Eigen::Vector3d &p_i_w_,
			Eigen::Vector3d &v_i_w_,
			Eigen::Quaterniond &q_i_w_,
			Eigen::Vector3d &omega_i,
			Eigen::Vector3d &a_i)
{
	p_i_w_ = p_i_w;
	v_i_w_ = v_i_w;
	q_i_w_ = q_i_w.toQuat();
	omega_i = omega_m_old - b_omega;
	a_i = a_m_old - b_a;
}

Eigen::Matrix<double,31,1>
EstimatorPredictor::GetStateVector(void)
{
	Eigen::Matrix<double,31,1> X;
	X << p_i_w, v_i_w, q_i_w.toQuat().coeffs(), b_omega, b_a, lambda, p_c_i, q_c_i.toQuat().coeffs(), p_ikf_w, q_ikf_w.toQuat().coeffs();
	return X;
}

EstimatorFull::EstimatorFull()
{
	// Covariance
	Matrix<double,28,1> P_vec;
	P_vec <<
	0, 0, 0, // P_p_i_w
	0, 0, 0, // P_v_i_w
	0, 0, 0, // P_q_i_w
	0, 0, 0, // P_b_omega;
	0, 0, 0, // P_b_a
	0, // P_lambda
	0, 0, 0, // P_p_c_i
	0, 0, 0, // P_q_c_i
	0, 0, 0, // P_p_ikf_w
	0, 0, 0; // P_q_ikf_w
	P = P_vec.asDiagonal();

	p_kf_v.setZero();
	q_kf_v.q.setIdentity();
}

void
EstimatorFull::SetCovarianceDiagonal(const Eigen::Matrix<double,28,1> &P_)
{
	P = P_.asDiagonal();
}

Eigen::Matrix<double,28,1>
EstimatorFull::GetCovarianceDiagonal(void)
{
	return P.diagonal();
}


void
EstimatorFull::PropagateCovariance(const Eigen::Vector3d &omega_m, const Eigen::Vector3d &a_m)
{
	//
	// Calculate state transition matrix
	//
	// Intermediate calculations
	double Delta_t2 = Delta_t*Delta_t/2.0;
	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();
	Matrix3d C_q_i_w_old = q_i_w_old.toQuat().toRotationMatrix();
	Matrix3d R = (C_q_i_w + C_q_i_w_old)/2.0;
	QuaternionAd q_i1_i( q_i_w.toQuat() * q_i_w_old.toQuat().conjugate() );
	Matrix3d R_i1_i = q_i1_i.toQuat().toRotationMatrix();
	Vector3d a_hat = (a_m - b_a);
	Vector3d a_hat_old = (a_m_old - b_a);
	Vector3d a_w_hat = C_q_i_w.transpose()*a_hat - g;
	Vector3d a_w_hat_old = C_q_i_w_old.transpose()*a_hat_old - g;
	Vector3d s = Delta_t/2.0 * ( R_i1_i.transpose()*a_hat + a_hat_old );
	Vector3d y = Delta_t2/2.0 * s;
	Matrix3d eye = Matrix3d::Identity();
	Matrix3d zero = Matrix3d::Zero();
	

	// State is:
	// x << p_i_w, v_i_w, q_i_w, b_omega, b_a, p_c_i, q_c_i, p_ikf_w, q_ikf_w
	// State transition matrix is:
	// F << F_p_p,  F_p_v,  F_p_q,  F_p_bomega,  F_p_ba,
	//      F_v_p,  F_v_v,  F_v_q,  F_v_bomega,  F_v_ba,
	//      F_q_p,  F_q_v,  F_q_q,  F_q_bomega,  F_q_ba,
	//      F_bo_p, F_bo_v, F_bo_q, F_bo_bomega, F_bo_ba,
	//      F_ba_p, F_ba_v, F_ba_q, F_ba_bomega, F_ba_ba,

	Matrix3d F_p_q, F_v_q, F_v_bo, F_p_bo, F_q_bo, F_p_ba, F_v_ba, F_q_ba;
	Matrix<double,15,15> F;

	F_p_q << -crossMat( C_q_i_w_old.transpose() * y );
	F_v_q << -crossMat( C_q_i_w_old.transpose() * s );

	F_v_bo << Delta_t2*(crossMat(a_w_hat))*R.transpose();
	F_p_bo << Delta_t/2.0 * F_v_bo;
	F_q_bo << -Delta_t*R.transpose();

	F_v_ba << -Delta_t*R.transpose();
	F_p_ba = Delta_t/2.0 * F_v_ba;
	F_q_ba << Matrix3d::Zero(); // Non zero if estimating gyro effect on accelerometer

	F << eye , eye*Delta_t, F_p_q, F_p_bo, F_p_ba,
	     zero, eye        , F_v_q, F_v_bo, F_v_ba,
	     zero, zero       , eye  , F_q_bo, F_q_ba,
	     zero, zero       , zero , eye   , zero  ,
	     zero, zero       , zero , zero  , eye   ;

	//
	// Calculate continous time noise matrix
	//
	Matrix<double,15, 15> N_c;
	N_c.setZero();
	for (int i = 0; i<3; i++)
	{
		N_c(3+i,3+i) = sq_sigma_a;
		N_c(6+i,6+i) = sq_sigma_omega;
		N_c(9+i,9+i) = sq_sigma_b_omega;
		N_c(12+i,12+i) = sq_sigma_b_a;
	}
	//
	// Calcuate Discrete noise matrix
	//
	Matrix<double,15,15> Q_d = Delta_t/2.0 * F * N_c * F.transpose() + Delta_t/2.0*N_c;

	//
	// Do the update
	//
	P.block<15,15>(0,0) = F * P.block<15,15>(0,0) * F.transpose() + Q_d;
	P.block(0,15,15,P.cols()-15) = F * P.block(0,15,15,P.cols()-15);
	P.block(15,0,P.rows()-15,15) = P.block(0,15,15,P.cols()-15).transpose();
}

void
EstimatorFull::UpdateCamera(const Eigen::Vector3d &p_c_v, const Eigen::Quaterniond &q_c_v,
	const Eigen::Matrix<double,6,6> &R,
	bool absolute, bool isKeyframe, double dLambda )
{
	static int i = 0;
	Eigen::QuaternionAd q_c_kf;
	Eigen::Vector3d p_c_kf;
	if (absolute)
	{
		q_c_kf = Eigen::QuaternionAd(q_c_v.conjugate() * q_kf_v.toQuat().conjugate());
	    p_c_kf = q_kf_v.toQuat().toRotationMatrix() * (p_c_v - p_kf_v);
	}
	else
	{
		q_c_kf.q = q_c_v;
		p_c_kf = p_c_v;
	}

	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();
	Matrix3d C_q_ikf_w = q_ikf_w.toQuat().toRotationMatrix();
	Matrix3d C_q_c_i = q_c_i.toQuat().toRotationMatrix();
	Matrix3d eye = Matrix3d::Identity();
	Matrix3d zero = Matrix3d::Zero();

	Matrix3d Hq_q_c_i, Hq_q_i_w, Hq_q_ikf_w;
	Matrix3d Hp_q_c_i, Hp_q_i_w, Hp_q_ikf_w, Hp_p_i_w, Hp_p_ikf_w, Hp_p_c_i;
	Vector3d Hp_lambda;
	Matrix<double,6,28> H;

	// Calculate measurement Jacobian
	Hp_p_i_w = C_q_c_i*C_q_ikf_w*exp(lambda);
	Hp_q_i_w = C_q_c_i*C_q_ikf_w*crossMat(C_q_i_w.transpose()*p_c_i)*exp(lambda);
	Hp_p_c_i = zero;C_q_c_i*C_q_ikf_w*C_q_i_w.transpose()*exp(lambda) - C_q_c_i*exp(lambda);
	Hp_q_c_i = zero;-C_q_c_i*crossMat(C_q_ikf_w*(p_i_w + C_q_i_w.transpose()*p_c_i - p_ikf_w) - p_c_i)*exp(lambda);
	Hp_p_ikf_w = -C_q_c_i*C_q_ikf_w*exp(lambda);
	Hp_q_ikf_w = -C_q_c_i*C_q_ikf_w*crossMat(p_i_w + C_q_i_w.transpose()*p_c_i - p_ikf_w)*exp(lambda);
	Hp_lambda = C_q_c_i*( C_q_ikf_w*( p_i_w + C_q_i_w.transpose()*p_c_i - p_ikf_w ) - p_c_i )*exp(lambda);

	Hq_q_i_w = C_q_c_i*C_q_i_w;
	Hq_q_c_i = C_q_c_i*(eye - C_q_i_w*C_q_ikf_w.transpose());
	Hq_q_ikf_w = -C_q_c_i*C_q_i_w;

	 H << Hp_p_i_w, zero, Hp_q_i_w, zero, zero,        Hp_lambda, Hp_p_c_i, Hp_q_c_i, Hp_p_ikf_w, Hp_q_ikf_w,
	          zero, zero, Hq_q_i_w, zero, zero, Vector3d::Zero(),     zero, Hq_q_c_i,       zero, Hq_q_ikf_w;

	// Calculate residual
	Matrix<double,6,1> r;
	Vector3d r_p = p_c_kf - C_q_c_i*( C_q_ikf_w*( p_i_w + C_q_i_w.transpose()*p_c_i - p_ikf_w ) - p_c_i )*exp(lambda);
	Quaterniond r_q_temp = q_c_kf.toQuat() * ( q_c_i.toQuat() * q_i_w.toQuat() * q_ikf_w.toQuat().conjugate() * q_c_i.toQuat().conjugate() ).conjugate();
	Vector3d r_q(-r_q_temp.x()*2.0,-r_q_temp.y()*2.0,-r_q_temp.z()*2.0);
	r_q = r_q;
	r << r_p, r_q;

	// Calculate kalman gain
	Matrix<double,6,6> S = H*P*H.transpose() + R;
	Matrix<double,28,6> K = P*H.transpose()*S.inverse();// (S.ldlt().solve(H*P)).transpose(); // P and S are symmetric
	Matrix<double,28,1>	x_error = K*r;

	// Apply Kalman gain
	i++;
	if ( i%1 == 0 ) 
	{
		/*std::clog << "P: " << P << std::endl <<
		"H: " << H << std::endl <<
		"K: " << K << std::endl <<
		"r: " << r.transpose() << std::endl;
*/
		ApplyCorrection( x_error );
		Matrix<double,28,28> T = Matrix<double,28,28>::Identity() - K*H;
		P = T*P*T.transpose() + K*R*K.transpose();
	}

	// Symmetrize
	P = (P + P.transpose())/2.0;

	if (isKeyframe)
	{
		// Add new keyframe from current position/camera pose
		UpdateKeyframe( );

		// Add noise to lambda
		P(15,15) += dLambda;

		// Update keyframe reference, if given path is absolute
		if (absolute)
		{
			p_kf_v = p_c_v;
			q_kf_v.q = q_c_v;
		}
	}
}

void inline 
EstimatorFull::ApplyCorrection( const Eigen::Matrix<double,28,1> &x_error ) 
{
	// TODO: look at quaternion error propagation. it smells
	p_i_w += x_error.segment<3>(0);
	v_i_w += x_error.segment<3>(3);	
	q_i_w = QuaternionAd(q_i_w.toQuat()*Quaterniond( 1, x_error(6)/2.0, x_error(7)/2.0, x_error(8)/2.0 ).conjugate()); q_i_w.normalize();
	b_omega += x_error.segment<3>(9);
	b_a += x_error.segment<3>(12);
	lambda += x_error(15);
	p_c_i += x_error.segment<3>(16);
	q_c_i = QuaternionAd(q_c_i.toQuat()*Quaterniond( 1, x_error(19)/2.0, x_error(20)/2.0, x_error(21)/2.0 ).conjugate()); q_c_i.normalize();
	p_ikf_w += x_error.segment<3>(22);
	q_ikf_w = QuaternionAd(q_ikf_w.toQuat()*Quaterniond( 1, x_error(25)/2.0, x_error(26)/2.0, x_error(27)/2.0 ).conjugate()); q_ikf_w.normalize();
}

void
EstimatorFull::UpdateKeyframe(void)
{
	Matrix<double,3,22> J_p;
	Matrix<double,3,22> J_q;
	Matrix<double,6,22> J;

	// Build jacobians
	J_p <<
	Matrix3d::Identity(),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Vector3d::Zero(),
	Matrix3d::Zero(),
	Matrix3d::Zero();

	J_q <<
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Matrix3d::Identity(),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Vector3d::Zero(),
	Matrix3d::Zero(),
	Matrix3d::Zero();

	J << J_p, J_q;

	// Calculate new keyframe pose
	p_ikf_w = p_i_w;
	q_ikf_w = q_i_w;

	// Calculate new P matrix
	P.block<22,6>(0,22) = P.block<22,22>(0,0)*J.transpose();
	P.block<6,6>(22,22) = J*P.block<22,22>(0,0)*J.transpose();
	P.block<6,22>(22,0) = P.block<22,6>(0,22).transpose();
}

}