﻿BatteryCapacityPredict
======================
利用多模型粒子滤波IMMPF对电池剩余容量进行跟踪和预测。（immpf_predict.m是主程序)
====================================================
2014/10/11
====================================================
添加了对失效时间的预测，但是由于多模型中的模型2的问题
无法在预定的时间内出现小于阈值的粒子，故模型2以及最后
的融合结果并没有实现。

2014/11/2    17:48 
====================================================
添加了权值的估计结果，并且融合了三个模型的估计

