#ifndef __LOCATION__
#define __LOCATION__

#define pi 3.1415926

Mat* filter(int Raw_data[], int length);
Mat* cal_center(Mat* pos);
Mat* delt_x(Mat* position);
void tick2angle(Mat* time);
Mat* ray(double angle1, double angle2);
Mat* sensor_distance(Mat* pos);
Mat* cal_mean(Mat* angles);
Mat* cal_init_mean(Mat* init_p);
Mat* Func(Mat* kk, Mat* iput);
Mat* Deriv(Mat* abc, Mat* iput, int n);
Mat* solve_x(Mat* h, Mat* y);
Mat* Func_R(Mat* kk, Mat* target_p);
Mat* Deriv_R(Mat* abc, Mat* iput, int n);
Mat* solve_R(Mat* real_p, Mat* target_p);
Mat* Func_r(Mat* kk, Mat* p1, Mat* p2);
Mat* Deriv_r(Mat* abc, Mat* p1, Mat* p2, int n);
Mat* solve_r(Mat* p1, Mat* p2);
Mat* Func_rotation_angle(Mat* angle_kk);
Mat* Deriv_rotation_angle(Mat* abc, int n);
Mat* solve_rotation_angle(Mat* pos, Mat* R, Mat* position, Mat* init);
int Outlier_detection(double* point);
Mat* cal_position(Mat* angle, Mat* pos);

#endif







