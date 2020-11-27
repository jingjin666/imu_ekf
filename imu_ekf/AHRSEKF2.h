#ifndef AHRSEKF2_H
#define AHRSEKF2_H

#include <vector>
#include "SensorData.h"

class EulerAngle
{


public:
	double Yaw;
	double Pitch;
	double Roll;

	EulerAngle(double yaw, double pitch, double roll)
	{
		Yaw = yaw;
		Pitch = pitch;
		Roll = roll;
	}
	EulerAngle()
	{
		Yaw = 0;
		Pitch = 0;
		Roll = 0;
	}
};
class AHRSEKF2
{
public:

	const double PI = 3.14159265358979323846;
	const double DEG_RAD = PI / 180;
	const double RAD_DEG = 180 / PI;

	AHRSEKF2();
	~AHRSEKF2();
	void ReadSensorData();
	void initalizevarMatrix(Eigen::Matrix<double, 7, 7> &PPrior0);
	void readSensorData();
	SensorData GetSensordatabyID(const long unsigned int &nId, bool flagnorm);
	EulerAngle InitializeEuler(const SensorData &sensordata);
	void InitializeVarMatrix(Eigen::Matrix<double, 7, 7> &Q, Eigen::Matrix<double, 3, 3> &R);
	void UpdateState(Eigen::Matrix<double, 1, 7> &x, Eigen::Matrix<double, 1, 7> &x_, const SensorData sensordata, const double T);
	void FillObserveState(Eigen::Matrix<double, 1, 3> &z, const SensorData sensordata);
	void FillTransiteMatrix(Eigen::Matrix<double, 7, 7> &Ak, const SensorData sensordata, Eigen::Matrix<double, 1, 7> &x, const double T);
	void FillObserveMatrix(const Eigen::Matrix<double, 1, 7> &x_, Eigen::Matrix<double, 1, 3> &hk, Eigen::Matrix<double, 3, 7> &Hk, const SensorData sensordata);
private:
	std::vector<SensorData> vSensorData;
};

#endif
