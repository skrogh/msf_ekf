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
	p_w_v = Vector3d(0, 0, 0);
	q_w_v = QuaternionAd(1, 0, 0, 0);

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
		const Eigen::Vector3d &p_w_v_,
		const Eigen::Quaterniond &q_w_v_)
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
	p_w_v = p_w_v_;
	q_w_v = Eigen::QuaternionAd( q_w_v_ );
	q_w_v.normalize();
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
	Vector3d v_i_w_old = v_i_w;
	Vector3d p_i_w_old = p_i_w;
	QuaternionAd q_i_w_old = q_i_w;

	// (eye(4) - Delta_t/4*Omega(omega_m - b_omega))*q_i_w_new == (eye(4) + Delta_t/4*Omega( omega_m_old - b_omega))*q_i_w_old
	Vector3d tmp = Delta_t/4.0*(omega_m - b_omega);
	Matrix4d Q1; // eye(4) - Delta_t/4*Omega(omega_m - b_omega)
	Q1 <<
	1, tmp(0), tmp(1), -tmp(2),
	-tmp(0), 1, tmp(2), tmp(1),
	-tmp(1), -tmp(2), 1, -tmp(0),
	tmp(2), -tmp(1), tmp(2), 1;

	tmp = Delta_t/4.0*(omega_m_old - b_omega);
	Matrix4d Q2; // eye(4) + Delta_t/4*Omega( omega_m_old - b_omega)
	Q2 <<
	1, tmp(0), tmp(1), -tmp(2),
	-tmp(0), 1, tmp(2), tmp(1),
	-tmp(1), -tmp(2), 1, -tmp(0),
	tmp(2), -tmp(1), tmp(2), 1;

	//q_i_w.coeffs() = Q1.lu().solve(Q2*q_i_w_old.coeffs()); 
	//q_i_w.normalize(); 
	tmp = Delta_t*(omega_m_old - b_omega)/2.0;
	Quaterniond Omega( 0, tmp(0), tmp(1), tmp(2) );
	q_i_w.q.coeffs() += (Omega.conjugate()*q_i_w.toQuat()).conjugate().coeffs();
	q_i_w.normalize();

	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();

	v_i_w += ( C_q_i_w_old.transpose()*(a_m_old - b_a) - g + C_q_i_w.transpose()*(a_m - b_a) - g )*Delta_t/2.0;
	p_i_w += (v_i_w_old + v_i_w)*Delta_t/2.0;

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
	X << p_i_w, v_i_w, q_i_w.toQuat().coeffs(), b_omega, b_a, lambda, p_c_i, q_c_i.toQuat().coeffs(), p_w_v, q_w_v.toQuat().coeffs();
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
	0, 0, 0, // P_p_w_v
	0, 0, 0; // P_q_w_v
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
	Vector3d omega_hat = omega_m - b_omega;
	Vector3d a_hat = a_m - b_a;

	Matrix3d C_q_w_v = q_w_v.toQuat().toRotationMatrix();
	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();
	Matrix3d C_q_c_i = q_c_i.toQuat().toRotationMatrix();

	double Delta_t2 = Delta_t*Delta_t/2.0;
	double Delta_t3 = Delta_t*Delta_t*Delta_t/6.0;
	double Delta_t4 = Delta_t*Delta_t*Delta_t*Delta_t/24.0;
	double Delta_t5 = Delta_t*Delta_t*Delta_t*Delta_t*Delta_t/120.0;

	// Build Jacobian
	Matrix3d A = -C_q_i_w.transpose()*crossMat(a_hat)
		*(Delta_t2*Matrix3d::Identity() - Delta_t3*crossMat(omega_hat) + Delta_t4*crossMat(omega_hat)*crossMat(omega_hat));
	Matrix3d B = -C_q_i_w.transpose()*crossMat(a_hat)
    	*(-Delta_t3*Matrix3d::Identity() + Delta_t4*crossMat(omega_hat) - Delta_t5*crossMat(omega_hat)*crossMat(omega_hat));
	Matrix3d C = -C_q_i_w.transpose()*crossMat(a_hat)
    	*(Delta_t*Matrix3d::Identity() - Delta_t2*crossMat(omega_hat) + Delta_t3*crossMat(omega_hat)*crossMat(omega_hat));
    Matrix3d D = -A;
	Matrix3d E = Matrix3d::Identity() - Delta_t*crossMat(omega_hat) + Delta_t2*crossMat(omega_hat)*crossMat(omega_hat);
	Matrix3d F = -(Delta_t*Matrix3d::Identity() - Delta_t2*crossMat(omega_hat) + Delta_t3*crossMat(omega_hat)*crossMat(omega_hat));
	Matrix<double,9,15> F_d_nonzero;

	// Optimized:
	Matrix<double,9,9> F_d_A;
	Matrix<double,9,6> F_d_B;

	Matrix<double,9,9> N_c_E = Matrix<double,9,9>::Zero();
	Matrix<double,6,6> N_c_H = Matrix<double,6,6>::Zero();

	Matrix<double,9,9> Q_d_A;
	Matrix<double,9,6> Q_d_B;
	Matrix<double,6,6> Q_d_C;

	#define P_C P.block<9,9>(0,0)
	#define P_D P.block<9,6>(0,9)
	#define P_E P.block<9,13>(0,9+6)
	#define P_D_ P.block<6,9>(9,0)
	#define P_E_ P.block<13,9>(9+6,0)
	#define P_G P.block<6,6>(9,9)
	#define P_H P.block<6,13>(9,9+6)
	



	F_d_A <<  Matrix3d::Identity(), Delta_t*Matrix3d::Identity(), A,
    Matrix3d::Zero(), Matrix3d::Identity(), C,
    Matrix3d::Zero(), Matrix3d::Zero(), E;

    F_d_B << B, -C_q_i_w.transpose()*Delta_t2,
    D, -C_q_i_w.transpose()*Delta_t,
    F, Matrix3d::Zero();

    for (int i = 0; i<3; i++)
    {
		N_c_E(3+i,3+i) = sq_sigma_a;
		N_c_E(6+i,6+i) = sq_sigma_omega;
		N_c_H(i,i) = sq_sigma_b_omega;
		N_c_H(3+i,3+i) = sq_sigma_b_a;
	}

	Q_d_A = Delta_t/2.0 * ( F_d_A*N_c_E*F_d_A.transpose() + F_d_B*N_c_H*F_d_B.transpose() + N_c_E );
	Q_d_B = Delta_t/2.0 * ( F_d_B*N_c_H );
	Q_d_C = Delta_t * N_c_H;

	Matrix<double,15,15> Q_d;
	Q_d << Q_d_A, Q_d_B, Q_d_B.transpose(), Q_d_C;

	
	P_C = F_d_A*P_C*F_d_A.transpose() + F_d_B*P_D.transpose()*F_d_A.transpose() + F_d_A*P_D*F_d_B.transpose() + F_d_B*P_G*F_d_B.transpose() + Q_d_A;
	P_D = F_d_A*P_D + F_d_B*P_G + Q_d_B;
	P_E = F_d_A*P_E + F_d_B*P_H;
	P_G += Q_d_C;
	P_D_ = P_D.transpose();
	P_E_ = P_E.transpose();
}

void
EstimatorFull::UpdateCamera(const Eigen::Vector3d &p_c_v, const Eigen::Quaterniond &q_c_v,
	const Eigen::Matrix<double,6,6> &R,
	bool absolute, bool isKeyframe, double dLambda )
{
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

	Matrix3d C_q_w_v = q_w_v.toQuat().toRotationMatrix();
	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();
	Matrix3d C_q_c_i = q_c_i.toQuat().toRotationMatrix();

	// Build Jacobian
	Matrix<double,3,28> H_p;
	H_p <<
	C_q_w_v.transpose() * exp(lambda),
	Matrix3d::Zero(),
	-C_q_w_v.transpose()*C_q_i_w.transpose()*crossMat(p_c_i)*exp(lambda),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	(C_q_w_v.transpose()*(C_q_i_w.transpose()*p_c_i + p_i_w) + p_w_v) * exp(lambda),
	C_q_w_v.transpose()*C_q_i_w.transpose()*exp(lambda),
	Matrix3d::Zero(),
	Matrix3d::Identity()*exp(lambda),
	-C_q_w_v.transpose()*crossMat( p_i_w + C_q_i_w.transpose()*p_c_i )*exp(lambda);

	Matrix<double,3,28> H_q;
	H_q <<
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	C_q_c_i,
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Vector3d::Zero(),
	Matrix3d::Zero(),
	Matrix3d::Identity(),
	Matrix3d::Zero(),
	C_q_c_i*C_q_i_w;

	/*Matrix3d::Zero(),
	Matrix3d::Zero(),
	-C_q_c_i,
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Vector3d::Zero(),
	Matrix3d::Zero(),
	-Matrix3d::Identity(),
	Matrix3d::Zero(),
	-C_q_c_i*C_q_i_w;*/

	Matrix<double,6,28> H;
	H << H_p, H_q;

	// Calculate residual

	Vector3d r_p = p_c_kf - ( C_q_w_v.transpose()*(p_i_w + C_q_i_w.transpose()*p_c_i) + p_w_v ) * exp(lambda);
	Quaterniond r_q = (q_c_kf.toQuat()*( q_c_i.toQuat()*q_i_w.toQuat()*q_w_v.toQuat() ).conjugate()).conjugate();

	Matrix<double,6,1> r;
	r << r_p, 2*r_q.x(), 2*r_q.y(), 2*r_q.z();

	// Calculate kalmn gain
	Matrix<double,6,6> S = H*P*H.transpose() + R;
	//K = P*H.transpose()*S^-1;
	//TODO: use cholsky to solve K*S = P*H.transposed()? 
	// (K*S)' = (P*H')'
	// S'*K' = H*P'
	// K' = S.transpose.ldlt().solve(H*P.transpose)
	Matrix<double,28,6> K = (S.ldlt().solve(H*P)).transpose(); // P and S are symmetric
	//Matrix<double,28,6> K = P*H.transpose()*S.inverse();//S.lu().solve(H*P).transpose();
	
	Matrix<double,28,1>	x_error = K*r;

	// Apply Kalman gain
	ApplyCorrection( x_error );
	Matrix<double,28,28> T = Matrix<double,28,28>::Identity() - K*H;
	P = T*P*T.transpose() + K*R*K.transpose();

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
	q_i_w = QuaternionAd(Quaterniond( 1, x_error(6)/2.0, x_error(7)/2.0, x_error(8)/2.0 ).conjugate()*q_i_w.toQuat()); q_i_w.normalize();
	b_omega += x_error.segment<3>(9);
	b_a += x_error.segment<3>(12);
	lambda += x_error(15);
	p_c_i += x_error.segment<3>(16);
	q_c_i = QuaternionAd(Quaterniond( 1, x_error(19)/2.0, x_error(20)/2.0, x_error(21)/2.0 ).conjugate()*q_c_i.toQuat()); q_c_i.normalize();
	p_w_v += x_error.segment<3>(22);
	q_w_v = QuaternionAd(Quaterniond( 1, x_error(25)/2.0, x_error(26)/2.0, x_error(27)/2.0 ).conjugate()*q_w_v.toQuat()); q_w_v.normalize();
}

void
EstimatorFull::UpdateKeyframe(void)
{
	
	Matrix3d C_q_i_w = q_i_w.toQuat().toRotationMatrix();
	Matrix3d C_q_c_i = q_c_i.toQuat().toRotationMatrix();

	Matrix<double,3,22> J_p;
	Matrix<double,3,22> J_q;
	Matrix<double,6,22> J;

	// Build jacobians
	J_p <<
	-C_q_c_i*C_q_i_w,
	Matrix3d::Zero(),
	C_q_c_i*crossMat( C_q_i_w*p_i_w ),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Vector3d::Zero(),
	-C_q_c_i,
	crossMat( C_q_c_i*(p_c_i + C_q_i_w*p_i_w) );

	J_q <<
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	-C_q_i_w.transpose(),
	Matrix3d::Zero(),
	Matrix3d::Zero(),
	Vector3d::Zero(),
	Matrix3d::Zero(),
	-C_q_c_i.transpose()*C_q_i_w.transpose();

	J << J_p, J_q;

	// Calculate new keyframe pose
	p_w_v = -C_q_c_i*( p_c_i + C_q_i_w*p_i_w );
	q_w_v = Eigen::QuaternionAd( q_i_w.toQuat().conjugate() * q_c_i.toQuat().conjugate() );

	// Calculate new P matrix
	P.block<22,6>(0,22) = P.block<22,22>(0,0)*J.transpose();
	P.block<6,6>(22,22) = J*P.block<22,22>(0,0)*J.transpose();
	P.block<6,22>(22,0) = P.block<22,6>(0,22).transpose();
}

}