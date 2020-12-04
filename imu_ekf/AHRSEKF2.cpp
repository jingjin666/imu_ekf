#include "AHRSEKF2.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Eigen/core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <math.h>

#include "Sensordata.h"
#include "Converter.h"

AHRSEKF2::AHRSEKF2()
{
}


AHRSEKF2::~AHRSEKF2()
{
}

void AHRSEKF2::initalizevarMatrix(Eigen::Matrix<double, 7, 7> &PPrior0)
{
	PPrior0 << 1, 0, 0, 0, 0, 0, 0,
					0, 1, 0, 0, 0, 0, 0,
					0, 0, 1, 0, 0, 0, 0,
					0, 0, 0, 1, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0,
					0, 0, 0, 0, 0, 0, 0;
}

void AHRSEKF2::readSensorData()
{
	printf("read the sensor raw data\n");
	SensorData sensordata;

	std::ifstream fin("session.txt");
	std::string line;
	int index = 0;
	while (std::getline(fin, line)) {
		std::istringstream sin(line);
		std::vector<std::string> fields;
		std::string field;
		while (getline(sin, field, ','))
		{
			fields.push_back(field);
		}

		sensordata.nId = index++;

		sensordata.Gyro.X = atof(fields[0].c_str());
		sensordata.Gyro.Y = atof(fields[1].c_str());
		sensordata.Gyro.Z = atof(fields[2].c_str());

		sensordata.Acc.X = atof(fields[3].c_str());
		sensordata.Acc.Y = atof(fields[4].c_str());
		sensordata.Acc.Z = atof(fields[5].c_str());

		sensordata.Mag.X = atof(fields[6].c_str());
		sensordata.Mag.Y = atof(fields[7].c_str());
		sensordata.Mag.Z = atof(fields[8].c_str());

		SensorData temp = sensordata;

		double norm = std::sqrt(temp.Acc.X*temp.Acc.X + temp.Acc.Y*temp.Acc.Y + temp.Acc.Z*temp.Acc.Z);
		temp.Acc.X /= norm;
		temp.Acc.Y /= norm;
		temp.Acc.Z /= norm;
		double pitch = asinf(temp.Acc.X);
		double roll = atan2f(temp.Acc.Y, temp.Acc.Z);

		double r1 = -sensordata.Mag.Y*cosf(roll) + sensordata.Mag.Z*sinf(roll);
		double r2 = sensordata.Mag.X*cosf(pitch) + sensordata.Mag.Y*sinf(pitch)*sinf(roll) + sensordata.Mag.Z*sinf(pitch)*cosf(roll);

		sensordata.EulerGroundTruth.Roll = roll;
		sensordata.EulerGroundTruth.Pitch = pitch;
		sensordata.EulerGroundTruth.Yaw = atan2(r1, r2) - 8.3 * DEG_RAD;

		vSensorData.push_back(sensordata);
	}

	fin.close();
}

void AHRSEKF2::ReadSensorData()
{
	std::cout << "read the sensor raw data" << std::endl;

	const unsigned long int ROW = 36, VOL = 4000;
	double d[VOL][ROW];
	std::ifstream in("myfile.txt");
	for (unsigned long int i = 0; i < VOL; i++)
	{
		for (int j = 0; j < ROW; j++)
		{
			in >> d[i][j];
		}
	}
	in.close();

	SensorData sensordata;
	for (unsigned long int i = 0; i < VOL; i++)
	{
		sensordata.nId = i;

		sensordata.Acc.X = d[i][8];
		sensordata.Acc.Y = d[i][9];
		sensordata.Acc.Z = d[i][10];

		sensordata.Gyro.X = d[i][26];
		sensordata.Gyro.Y = d[i][27];
		sensordata.Gyro.Z = d[i][28];

		sensordata.Mag.X = d[i][14];
		sensordata.Mag.Y = d[i][15];
		sensordata.Mag.Z = d[i][16];

		sensordata.EulerGroundTruth.Roll = d[i][29];
		sensordata.EulerGroundTruth.Pitch = d[i][30];
		sensordata.EulerGroundTruth.Yaw = d[i][31];

		vSensorData.push_back(sensordata);
	}

	//std::cout.precision(10);
	//for (int i = 0; i < 10;i++)
	//{
	//	std::cout << vSensorData.at(i).Gyro.X << std::endl;
	//}

	std::cout << "finish loading the dataset" << std::endl;
}

SensorData AHRSEKF2::GetSensordatabyID(const long unsigned int &nId, bool flagnorm)
{
	SensorData sensordata = vSensorData.at(nId);

	if (flagnorm == true)
	{
		double norm = std::sqrt(sensordata.Acc.X*sensordata.Acc.X + sensordata.Acc.Y*sensordata.Acc.Y + sensordata.Acc.Z*sensordata.Acc.Z);
		sensordata.Acc.X /= norm;
		sensordata.Acc.Y /= norm;
		sensordata.Acc.Z /= norm;

		norm = std::sqrt(sensordata.Gyro.X*sensordata.Gyro.X + sensordata.Gyro.Y*sensordata.Gyro.Y + sensordata.Gyro.Z*sensordata.Gyro.Z);
		sensordata.Gyro.X /= norm;
		sensordata.Gyro.Y /= norm;
		sensordata.Gyro.Z /= norm;

		norm = std::sqrt(sensordata.Mag.X*sensordata.Mag.X + sensordata.Mag.Y*sensordata.Mag.Y + sensordata.Mag.Z*sensordata.Mag.Z);
		sensordata.Mag.X /= norm;
		sensordata.Mag.Y /= norm;
		sensordata.Mag.Z /= norm;
	}
	else;

	return sensordata;
}

EulerAngle AHRSEKF2::InitializeEuler(const SensorData &sensordata)
{
	double pitch, roll, yaw;
	SensorData temp = sensordata;

	double norm = std::sqrt(temp.Acc.X*temp.Acc.X + temp.Acc.Y*temp.Acc.Y + temp.Acc.Z*temp.Acc.Z);
	temp.Acc.X /= norm;
	temp.Acc.Y /= norm;
	temp.Acc.Z /= norm;
	pitch = asinf(temp.Acc.X);
	roll = atan2f(temp.Acc.Y, temp.Acc.Z);

	double r1 = -temp.Mag.Y*cos(roll) + temp.Mag.Z*sin(roll);
	double r2 = temp.Mag.X*cos(pitch) + temp.Mag.Y*sin(pitch)*sin(roll) + temp.Mag.Z*sin(pitch)*cos(roll);

	yaw = atan2(r1, r2) - 8.3 * DEG_RAD;

	return EulerAngle(yaw, pitch, roll);
}

// Q 7*7(4*4 1e-6 单位阵 3*3 1e-8 单位阵) R 6*6(0.1 单位阵)
void AHRSEKF2::InitializeVarMatrix(Eigen::Matrix<double, 7, 7> &Q, Eigen::Matrix<double, 3, 3> &R)
{
    //过程噪声
	const double w_process_noise_var = 0.0001;  // rot vel var 
	const double wb_process_noise_var = 0.003; // gyro bias change var

    //测量噪声
	const double acc_measure_noise_var = 0.03;  // acc var 
	
	Q = Eigen::MatrixXd::Identity(7, 7);
	Q.block<4, 4>(0, 0) *= w_process_noise_var;
	Q.block<3, 3>(4, 4) *= wb_process_noise_var;

	R = Eigen::MatrixXd::Identity(3, 3);
	R.block<3, 3>(0, 0) *= acc_measure_noise_var;

	//std::cout << Q << std::endl;
	//std::cout << R << std::endl;
}

// x q.w q.x q.y q.z bx by bz
//   0   1   2   3   4  5  6
void AHRSEKF2::UpdateState(Eigen::Matrix<double, 1, 7> &x, Eigen::Matrix<double, 1, 7> &x_,const SensorData sensordata, const double T)
{
	double Gyro_Xcorrect, Gyro_Ycorrect, Gyro_Zcorrect;
	
	Gyro_Xcorrect = sensordata.Gyro.X - x[4];
	Gyro_Ycorrect = sensordata.Gyro.Y - x[5];
	Gyro_Zcorrect = sensordata.Gyro.Z - x[6];

	x_[0] = x[0] + (-x[1]*Gyro_Xcorrect - x[2]*Gyro_Ycorrect - x[3]*Gyro_Zcorrect) * T/2;
	x_[1] = x[1] + ( x[0]*Gyro_Xcorrect + x[2]*Gyro_Zcorrect - x[3]*Gyro_Ycorrect) * T/2;
	x_[2] = x[2] + ( x[0]*Gyro_Ycorrect - x[1]*Gyro_Zcorrect + x[3]*Gyro_Xcorrect) * T/2;
	x_[3] = x[3] + ( x[0]*Gyro_Zcorrect + x[1]*Gyro_Ycorrect - x[2]*Gyro_Xcorrect) * T/2;

	//printf("%f, %f, %f, %f\n", x_[0], x_[1], x_[2], x_[3]);

	x_[4] = x[4];
	x_[5] = x[5];
	x_[6] = x[6];

	double norm;
	norm = sqrt(x_[0]*x_[0] + x_[1]*x_[1] + x_[2]*x_[2] + x_[3]*x_[3]);
	x_[0] /= norm;
	x_[1] /= norm;
	x_[2] /= norm;
	x_[3] /= norm;

	//printf("%f, %f, %f, %f\n", x_[0], x_[1], x_[2], x_[3]);
}

void AHRSEKF2::FillObserveState(Eigen::Matrix<double, 1, 3> &z, const SensorData sensordata)
{
	z[0] = sensordata.Acc.X;
	z[1] = sensordata.Acc.Y;
	z[2] = sensordata.Acc.Z;
}

void AHRSEKF2::FillTransiteMatrix(Eigen::Matrix<double, 7, 7> &Ak, const SensorData sensordata, Eigen::Matrix<double, 1, 7> &x, const double T)
{
	Ak = Eigen::MatrixXd::Zero(7, 7);

	Ak.block<4, 4>(0, 0) << 0, -(sensordata.Gyro.X - x[4]), -(sensordata.Gyro.Y - x[5]), -(sensordata.Gyro.Z - x[6]),
							(sensordata.Gyro.X - x[4]), 0, (sensordata.Gyro.Z - x[6]), -(sensordata.Gyro.Y - x[5]),
							(sensordata.Gyro.Y - x[5]), -(sensordata.Gyro.Z - x[6]), 0, (sensordata.Gyro.X - x[4]),
							(sensordata.Gyro.Z - x[6]), (sensordata.Gyro.Y - x[5]), -(sensordata.Gyro.X - x[4]), 0;
	
	Ak.block<4, 3>(0, 4) <<  x[1],  x[2],  x[3],
							-x[0],  x[3], -x[2],
							-x[3], -x[0],  x[1],
							 x[2], -x[1], -x[0];

	Ak = Eigen::MatrixXd::Identity(7, 7) + 0.5 * T * Ak;
}

void  AHRSEKF2::FillObserveMatrix(const Eigen::Matrix<double, 1, 7> &x_, Eigen::Matrix<double, 1, 3> &hk, Eigen::Matrix<double, 3, 7> &Hk, const SensorData sensordata)
{
	Hk.block<3, 7>(0, 0) << -x_[2], x_[3], -x_[0], x_[1], 0, 0, 0,
							x_[1],  x_[0],  x_[3], x_[2], 0, 0, 0,
							x_[0], -x_[1], -x_[2], x_[3], 0, 0, 0;

	Hk = 2 * Hk;

	hk = x_ * Hk.transpose();
}

// Hk
//void AHRSEKF2::FillStateGain(Eigen::Matrix<double, 7, 7>)
//{
//
//}