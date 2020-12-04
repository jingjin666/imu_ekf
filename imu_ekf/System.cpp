#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
using namespace std;

#include "System.h"
#include "AHRSEKF2.h"
#include "SensorData.h"
#include "Converter.h"
#include <stdio.h>

System::System()
{
}


System::~System()
{
}

int System::RunEKF2()
{
	AHRSEKF2 ekf;
	const double T = 0.01;
	unsigned int index = 0;
	SensorData sensordatanormk;
	SensorData sensordataunnormk;

	ekf.readSensorData();

	Eigen::Matrix<double, 1, 7> x = Eigen::MatrixXd::Zero(1, 7);
	Eigen::Matrix<double, 1, 7> x_ = Eigen::MatrixXd::Zero(1, 7);
	Eigen::Matrix<double, 1, 3> z = Eigen::MatrixXd::Zero(1, 3);

	Eigen::Matrix<double, 7, 7> Pk_ = Eigen::MatrixXd::Zero(7, 7);
	Eigen::Matrix<double, 7, 7> Pk = Eigen::MatrixXd::Zero(7, 7);

	Eigen::Matrix<double, 1, 3> hk;
	Eigen::Matrix<double, 3, 7> Hk = Eigen::MatrixXd::Zero(3, 7);
	Eigen::Matrix<double, 7, 3> Kk = Eigen::MatrixXd::Zero(7, 3);

	Eigen::Matrix<double, 3, 3> R;
	Eigen::Matrix<double, 7, 7> Q;
	Eigen::Matrix<double, 7, 7> Ak;

	Eigen::Vector3d euler;
	Eigen::Quaterniond qfilter;

    //用加速度计数据初始化欧拉角
	EulerAngle eulerinit;
	eulerinit = ekf.InitializeEuler(ekf.GetSensordatabyID(0,false));
    //初始化欧拉角转换成四元数矩阵,使得系统更快收敛
	Eigen::Quaterniond qinit = Converter::euler2quat(Eigen::Vector3d(eulerinit.Yaw,eulerinit.Pitch,eulerinit.Roll));
    x[0] = qinit.w(), x[1] = qinit.x(), x[2] = qinit.y(), x[3] = qinit.z(), x[4] = 0.01, x[5] = 0.01, x[6] = 0.01;

	float roll_init = atan2f(2.0f*(x[0] * x[1] + x[2] * x[3]), 1.0f - 2.0f*(x[1] * x[1] + x[2] * x[2])) * ekf.RAD_DEG;
	float pitch_init = -asinf(2.0f*(x[1]* x[3] - x[0] * x[2])) * ekf.RAD_DEG;
	printf("init_roll is %f, %f\n", roll_init, pitch_init);

    //初始化协方差矩阵
	ekf.initalizevarMatrix(Pk);

    //初始化过程噪声和测量噪声矩阵
	ekf.InitializeVarMatrix(Q,R);
    
	//Eigen::Matrix<double, 7, 3> s = Eigen::MatrixXd::Zero(7, 3);
	//Eigen::Matrix<double, 7, 3> h_t = Eigen::MatrixXd::Zero(7, 3);
	//Eigen::Matrix<double, 3, 3> sp = Eigen::MatrixXd::Zero(3, 3);
	//Eigen::Matrix<double, 3, 3> sp_i = Eigen::MatrixXd::Zero(3, 3);
	//Eigen::Matrix<double, 1, 3> z_ = Eigen::MatrixXd::Zero(1, 3);
	//Eigen::Matrix<double, 3, 1> z_t = Eigen::MatrixXd::Zero(3, 1);
	
	ofstream fout("C:\\Users\\Administrator\\PycharmProjects\\untitled\\EKF_OUTPUT.TXT");
	while(1)
	{
		//  the sensordatak that have the normalized, if want to use the unnormalized data, set the false flag in the GetSensordata function  
		// 归一化的数据
		sensordatanormk = ekf.GetSensordatabyID(index,true);
        // 没有归一化的数据
		sensordataunnormk = ekf.GetSensordatabyID(index,false);

		//填充观测数据Z，归一化后的加速度计数据
		ekf.FillObserveState(z,sensordatanormk);
		
		//预测状态矩阵X
		ekf.UpdateState(x,x_,sensordataunnormk,T);
		//printf("Pre:%f, %f, %f, %f\n", x_[0], x_[1], x_[2], x_[3]);

		//获取状态转移矩阵F
		ekf.FillTransiteMatrix(Ak, sensordataunnormk, x, T);

		//预测协方差矩阵P
		Pk_ = Ak * Pk * Ak.transpose() + Q;

		//计算残差hk和观测矩阵Hk
		ekf.FillObserveMatrix(x_,hk,Hk,sensordatanormk);

		//计算卡尔曼增益
		//h_t = Hk.transpose();
		//s = Pk_ * h_t;
		//sp = Hk * Pk_ * h_t + R;
		//sp_i = sp.inverse();
		//Kk = s * sp_i;
		Kk = Pk_ * Hk.transpose() * (Hk * Pk_ * Hk.transpose() + R).inverse();

        //更新状态矩阵
		//z_ = z - hk;
		//z_t = z_.transpose();
		//x = (x_.transpose() + Kk * z_t).transpose();
		x = (x_.transpose() + Kk * (z - hk).transpose()).transpose();
		//printf("Update:%f, %f, %f, %f\n", x[0], x[1], x[2], x[3]);

        //归一化四元数
		Converter::Normalize(x);
		//printf("Update:%f, %f, %f, %f\n", x[0], x[1], x[2], x[3]);

        //输出滤波后的欧拉角
		qfilter = Converter::vector4d2quat(x.block<1,4>(0, 0));
		euler = Converter::quat2euler(qfilter);// yaw pitch roll
		euler[0] = euler[0] * ekf.RAD_DEG  - 8.3;
		euler[1] = euler[1] * ekf.RAD_DEG;
		euler[2] = euler[2] * ekf.RAD_DEG;
        //std::cout << "++++euler:" << euler.transpose() << std::endl;
        
		//输出原始加速度计算的欧拉角
        //float acc_roll = -atan2f(sensordataunnormk.Acc.X, sqrt(sensordataunnormk.Acc.Y*sensordataunnormk.Acc.Y + sensordataunnormk.Acc.Z*sensordataunnormk.Acc.Z)) * ekf.RAD_DEG;
        //float acc_pitch = atan2f(sensordataunnormk.Acc.Y, sqrt(sensordataunnormk.Acc.X*sensordataunnormk.Acc.X + sensordataunnormk.Acc.Z*sensordataunnormk.Acc.Z)) * ekf.RAD_DEG;
		sensordatanormk.EulerGroundTruth.Yaw *= ekf.RAD_DEG;
		sensordatanormk.EulerGroundTruth.Pitch *= ekf.RAD_DEG;
		sensordatanormk.EulerGroundTruth.Roll *= ekf.RAD_DEG;
		//std::cout << "-----euler:"  << sensordatanormk.EulerGroundTruth.Roll << " " << sensordatanormk.EulerGroundTruth.Pitch << " "  << sensordatanormk.EulerGroundTruth.Yaw << std::endl; 

		//printf("%f,%f\n", euler[2], sensordatanormk.EulerGroundTruth.Roll);
		//fout << euler[2] << "," << sensordatanormk.EulerGroundTruth.Roll << endl;
		fout << euler[1] << "," << sensordatanormk.EulerGroundTruth.Pitch << endl;

        //更新协方差矩阵
		Pk = (Eigen::MatrixXd::Identity(7, 7) - Kk * Hk) * Pk_;

		index++;
		if (index == 3000) {
			break;
		}
	}
	fout.close();

	return 0;
}

int System::SaveData(std::vector<SensorData> vSensorData)
{
	std::ofstream outfile;
	long unsigned int index = 0;;
	outfile.open("log.txt");
	if (!outfile)
	{
		std::cout << "sava data fail" << std::endl;
		return 0;
	}

	while(index < (4200-2))
	{
		outfile << index;
		outfile << " ";
		outfile << vSensorData.at(index).CalculateEuler.Yaw;
		outfile << " ";
		outfile << vSensorData.at(index).CalculateEuler.Pitch;
		outfile << " ";
		outfile << vSensorData.at(index).CalculateEuler.Roll;
		outfile << " ";
		outfile << vSensorData.at(index).EulerGroundTruth.Yaw;
		outfile << " ";
		outfile << vSensorData.at(index).EulerGroundTruth.Pitch;
		outfile << " ";
		outfile << vSensorData.at(index).EulerGroundTruth.Roll;
		outfile << std::endl;

		index++;

	}

	outfile.close();

	return 0;
}
