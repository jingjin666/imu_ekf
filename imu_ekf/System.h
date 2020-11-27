#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include "SensorData.h"

class System
{
public:
	System();
	~System();
	int RunEKF();
	int RunEKF2();
	int RunESKF();
	int SaveData(std::vector<SensorData> vSensorData);
};
#endif
