#include <stdio.h>
#include<math.h>
#include<time.h>
#include "light_matrix.h"
#include "location.h"

// 数据预处理，提取有效数据，并返回矩阵
Mat* filter(int Raw_data[],int length){
    // 如果数组长度小于72，即有效收到的数据量小于3，则返回NULL
    if (length < 72)
        return NULL;
    // 计数，用于测试是否为第一组数据
    int count = 0;
    // 接受数据的缓存
    Mat filter_data;
    // 返回的有效数据
    static Mat result;
    MatCreate(&filter_data, 1, 8);
    for (int i = 0; i < length;) {
        if (Raw_data[i] != 0xf0)
            i++;
        else{
            // 检测对定位置的标志位是否对应
            if (Raw_data[i] == 0xf0 && Raw_data[i + 3] == 0xf2 && Raw_data[i + 6] == 0xf4 && Raw_data[i + 9] == 0xf6 &&
                Raw_data[i + 12] == 0xf1 && Raw_data[i + 15] == 0xf3 && Raw_data[i + 18] == 0xf5 && Raw_data[i + 21] == 0xf7) {
                for (int j=0; j < 8; j++)
                    filter_data.element[0][j] = 256.0 * Raw_data[i + 3 * j + 1] + Raw_data[i + 3 * j + 2];
                if (count == 0) {
                    MatCopy(&filter_data, &result);
                    count++;
                }
                else {
                    MatCopy(MatRowStock(&result, &filter_data),&result);
                }
                i += 24;
            }
            else
                i++;
        }
    }
    // 释放内存
    MatDelete(&filter_data);
    return &result;
}

// 计算中心坐标位置
Mat* cal_center(Mat* pos) {
    static Mat point;
    MatCreate(&point, 1, 3);
    for (int i = 0; i < 3; i++) {
        point.element[0][i] = MatSum(MatCol(pos, i)) / pos->row;
    }
    return &point;
}

// 去除平移量，只保留旋转量
Mat* delt_x(Mat* position){
    Mat center;
    static Mat delt_position;
    MatCreate(&center, 1, 3);
    MatCreate(&delt_position, position->row, position->col);
    MatCopy(cal_center(position), &center);
    for (int i=0;i<position->row;i++)
        for (int j = 0; j < position->col; j++) {
            delt_position.element[i][j] = position->element[i][j] - center.element[0][j];
        }
    MatDelete(&center);
    return &delt_position; 
}

// 时间转换为角度
void tick2angle(Mat* time){
    int row, col;
    for (row = 0; row < time->row; row++)
        for (col = 0; col < time->col; col++)
            time->element[row][col] = pi * (time->element[row][col] / 1.20 - 4000.0) / 8333.0;
}

// 返回横扫与竖扫平面相交的线向量（基站到传感器的单位线向量）
Mat* ray(double angle1,double angle2){
    Mat plane1, plane2;
    static Mat res;
    MatCreate(&plane1, 1, 3);
    MatCreate(&plane2, 1, 3);
    double val1[3] = { cos(angle1),0,-1 * sin(angle1) };
    double val2[3] = { 0,cos(angle2), sin(angle2) };
    MatSetVal(&plane1, val1);
    MatSetVal(&plane2, val2);
    MatCopy(MatUnit(MatCross3(&plane2, &plane1)), &res);
    MatDelete(&plane1);
    MatDelete(&plane2);
    return &res;

}

// 计算两个传感器之间的距离
Mat* sensor_distance(Mat* pos){
    static Mat dist;
    MatCreate(&dist, 1, 6);
    
    int count = 0;
    for(int row=0;row<pos->row;row++)
        for (int j = row+1; j < pos->row; j++) {
            dist.element[0][count] = sqrt(MatSum(MatPow(MatSub(MatRow(pos, row), MatRow(pos, j)), 2)));
            count += 1;
        }
    return &dist;
} 

// 计算角度的平均值
Mat* cal_mean(Mat* angles) {
    int  col;
    static Mat mean;
    MatCreate(&mean, 1, angles->col);
    for (col = 0; col < angles->col; col++) {
        mean.element[0][col] = MatSum(MatCol(angles, col)) / angles->row;
    }
    return &mean;
}

// 计算初始坐标的平均值
Mat* cal_init_mean(Mat* init_p) {
    static Mat init_mean;
    MatCreate(&init_mean, 4, 3);
    double t;
    t = init_p->element[1][1];
    for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
            init_mean.element[j][k] = 0.2 * ((init_p + 0)->element[j][k] + (init_p + 1)->element[j][k] +
                (init_p + 2)->element[j][k] + (init_p + 3)->element[j][k] +
                (init_p + 4)->element[j][k]);
    return &init_mean;
}
//LM最优化算法对方程进行求解
// 位置函数
Mat* Func(Mat* kk, Mat* iput) {
    Mat k;
    if (kk->row != 1) {
        MatCopy(MatTrans(kk), &k);
    }
    else {
        MatCopy(kk,&k);
    }
    static Mat result;
    MatCreate(&result, 1, 6);
    int count = 0;
    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++) {
            result.element[0][count] = k.element[0][i] * k.element[0][i] + k.element[0][j] * k.element[0][j] -
                2 * k.element[0][i] * k.element[0][j] * iput->element[0][count];
            count += 1;
        }
    MatDelete(&k);
    return &result;
}
// 位置函数微分
Mat* Deriv(Mat* abc, Mat* iput, int n) {
    Mat x1, x2;
    Mat p1,  p2;
    static Mat d;
    MatCopy(abc, &x1);
    MatCopy(abc, &x2);
    x1.element[n][0] -= 0.00001;
    x2.element[n][0] += 0.00001;
    MatCopy(Func(&x1, iput),&p1);
    MatCopy(Func(&x2, iput),&p2);
    MatCopy(MatTrans(MatNumMul(MatSub(&p2, &p1), 50000.0)),&d);
    MatDelete(&x1);
    MatDelete(&x2);
    MatDelete(&p1);
    MatDelete(&p2);
    return &d;
}
// 计算四个传感器与基站之间的距离
Mat* solve_x(Mat* h,Mat* y){
    // 设置迭代计数，总迭代数，方程个数
    int step = 0, conve = 40;
    double mse, mse_tmp,q;
    // 雅克比矩阵，误差矩阵
    Mat J,fx,fx_tmp,H,dx,xk_tmp;
    MatCreate(&J, 6, 4);
    MatCreate(&fx, 6, 1);
    MatCreate(&fx_tmp, 6, 1);
    MatCreate(&H, 3, 3);
    // 上一次迭代的均方差
    double lase_mse = 0.0;
    double u = 1, v = 2;
    //初始值设置
    Mat xk;
    static Mat k;
    MatCreate(&xk, 4, 1);
    MatCreate(&k, 1, 4);
    double buff[] = { 1500.0,1500.0,1500.0,1500.0 };
    MatSetVal(&xk, buff);
    MatSetVal(&k, buff);
    while (conve) {
        mse = 0.0;
        mse_tmp = 0.0;
        step += 1;
        MatCopy(MatSub(Func(&k, h), y),&fx);
        mse += MatSum(MatPow(&fx, 2));
        for (int j = 0; j < 4; j++) {
            MatColStead(&J, Deriv(&xk, h, j), j);
        }
        
        MatCopy(MatAdd(MatMul(MatTrans(&J), &J), MatNumMul(MatEye(4), u)),&H);
        MatCopy(MatMul(MatMul(MatInv(MatNumMul(&H, -1)), MatTrans(&J)), MatTrans(&fx)),&dx);  
        MatCopy(&xk, &xk_tmp);    
        MatCopy(MatAdd(&xk_tmp, &dx), &xk_tmp);
        MatCopy(MatSub(Func(&xk_tmp, h), y), &fx_tmp);
        mse_tmp = MatSum(MatPow(MatCol(&fx_tmp, 0), 2));
        // 判断是否下降
        q = (mse - mse_tmp) / (MatNumMul(MatMul(MatTrans(&dx) ,MatSub(MatNumMul(&dx,u) , MatMul(MatTrans(&J) , MatTrans(&fx)))),0.5)->element[0][0]);
        if (q > 0) {
            double s = 1.0 / 3;
            v = 2;
            mse = mse_tmp;
            MatCopy(&xk_tmp, &xk);
            MatCopy(MatTrans(&xk), &k);
            double temp = 1.0 - (2 * q - 1) * (2 * q - 1) * (2 * q - 1);
            if (s > temp)
                u = u * s;
            else
                u = u * temp;
        }
        else {
            u = u * v;
            v = 2 * v;
            MatCopy(&xk_tmp, &xk);
            MatCopy(MatTrans(&xk), &k);
        }
        if (mse - lase_mse<0.001 && mse - lase_mse>-0.001) {
            MatDelete(&J);
            MatDelete(&fx);
            MatDelete(&fx_tmp);
            MatDelete(&H);
            MatDelete(&dx);
            MatDelete(&xk_tmp);
            MatDelete(&xk);
            return &k;
        }
            
        //printf("error:%f\n", (mse - lase_mse));
        lase_mse = mse;
        conve -= 1;
    }
    return NULL;
}


//计算旋转矩阵1
Mat* Func_R(Mat* kk, Mat* target_p) {
    Mat k,trans_p;
    if (kk->row != 1) {
        MatCopy(MatTrans(kk), &k);
    }
    else {
        MatCopy(kk, &k);
    }
    MatCopy(MatTrans(target_p), &trans_p);
    static Mat result;
    MatCreate(&result, 1, 15);
    int count = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
        {
            result.element[0][count] = k.element[0][2 * i + 0] * trans_p.element[0][j] + k.element[0][2 * i + 1] * trans_p.element[1][j];
            count += 1;
        }
    result.element[0][count] = k.element[0][0] * k.element[0][0] + k.element[0][2] * k.element[0][2] + k.element[0][4] * k.element[0][4];
    result.element[0][count + 1] = k.element[0][1] * k.element[0][1] + k.element[0][3] * k.element[0][3] + k.element[0][5] * k.element[0][5];
    result.element[0][count + 2] = k.element[0][1] * k.element[0][0] + k.element[0][3] * k.element[0][2] + k.element[0][5] * k.element[0][4];
    MatDelete(&k);
    MatDelete(&trans_p);

    return &result;
}
// 微分运算
Mat* Deriv_R(Mat* abc, Mat* iput, int n) {
    Mat x1, x2;
    Mat p1, p2;
    static Mat d;
    MatCopy(abc, &x1);
    MatCopy(abc, &x2);
    x1.element[n][0] -= 0.000001;
    x2.element[n][0] += 0.000001;
    MatCopy(Func_R(&x1, iput),&p1);
    MatCopy(Func_R(&x2, iput),&p2);
    MatCopy(MatTrans(MatNumMul(MatSub(&p2,&p1), 500000.0)),&d);
    MatDelete(&x1);
    MatDelete(&x2);
    MatDelete(&p1);
    MatDelete(&p2);
    return &d;
}

Mat* solve_R(Mat* real_p, Mat* target_p) {
    Mat y,delt_real_p;
    MatCreate(&delt_real_p, real_p->row, real_p->col);
    MatCopy(delt_x(real_p), &delt_real_p);
    MatCreate(&y, 1, 15);
    int count = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++) {
            y.element[0][count] = delt_real_p.element[j][i];
            count += 1;
        }
    y.element[0][count] = 1;
    y.element[0][count + 1] = 1;
    y.element[0][count + 2] = 0;
    // 设置迭代计数，总迭代数，方程个数
    int step = 0, conve = 40;
    double mse, mse_tmp, q;
    // 雅克比矩阵，误差矩阵
    Mat J, fx, fx_tmp, H, dx, xk_tmp;
    MatCreate(&J, 15, 6);
    MatCreate(&fx, 15, 1);
    MatCreate(&fx_tmp, 15, 1);
    MatCreate(&H, 3, 3);
    // 上一次迭代的均方差
    double lase_mse = 0.0;
    double u = 1, v = 2;
    //初始值设置
    Mat xk;
    static Mat k;
    MatCreate(&xk, 6, 1);
    MatCreate(&k, 1, 6);
    double buff[] = { 1,1,1,1,1,1 };
    MatSetVal(&xk, buff);
    MatSetVal(&k, buff);
    while (conve) {
        mse = 0.0;
        mse_tmp = 0.0;
        step += 1;
        MatCopy(MatSub(Func_R(&k, target_p), &y), &fx);
        mse += MatSum(MatPow(&fx, 2));
        for (int j = 0; j < 6; j++) {
            MatColStead(&J, Deriv_R(&xk, target_p, j), j);
        }
        MatCopy(MatAdd(MatMul(MatTrans(&J), &J), MatNumMul(MatEye(6), u)), &H);
        MatCopy(MatMul(MatMul(MatInv(MatNumMul(&H, -1)), MatTrans(&J)), MatTrans(&fx)), &dx);
        MatCopy(&xk, &xk_tmp);
        MatCopy(MatAdd(&xk_tmp, &dx), &xk_tmp);
        MatCopy(MatSub(Func_R(&xk_tmp, target_p), &y), &fx_tmp);
        mse_tmp = MatSum(MatPow(MatCol(&fx_tmp, 0), 2));
        // 判断是否下降
        q = (mse - mse_tmp) / (MatNumMul(MatMul(MatTrans(&dx), MatSub(MatNumMul(&dx, u), MatMul(MatTrans(&J), MatTrans(&fx)))), 0.5)->element[0][0]);
        if (q > 0) {
            double s = 1.0 / 3;
            v = 2;
            mse = mse_tmp;
            MatCopy(&xk_tmp, &xk);
            MatCopy(MatTrans(&xk), &k);
            double temp = 1.0 - (2 * q - 1) * (2 * q - 1) * (2 * q - 1);
            if (s > temp)
                u = u * s;
            else
                u = u * temp;
        }
        else {
            u = u * v;
            v = 2 * v;
            MatCopy(&xk_tmp, &xk);
            MatCopy(MatTrans(&xk), &k);
        }
        if (mse - lase_mse<0.00001 && mse - lase_mse>-0.00001) {
            Mat  xy, z;
            static Mat result;
            MatCreate(&result, 3, 3);
            MatCreate(&xy, 3, 2);
            MatCreate(&z, 3, 1);
            MatSetVal(&xy, k.element[0]);
            MatCopy(MatNumMul(MatCross3(MatTrans(MatCol(&xy, 0)), MatTrans(MatCol(&xy, 1))), -1), &z);
            MatCopy(MatColStock(&xy, MatTrans(&z)),&result);

            MatDelete(&y);
            MatDelete(&k);
            MatDelete(&delt_real_p);
            MatDelete(&J);
            MatDelete(&fx);
            MatDelete(&fx_tmp);
            MatDelete(&H);
            MatDelete(&dx);
            MatDelete(&xk_tmp);
            MatDelete(&xk);
            MatDelete(&xy);
            MatDelete(&z);
            return &result;
        }
            
        //printf("error:%f\n", (mse - lase_mse));
        lase_mse = mse;
        conve -= 1;
    }
    return NULL;
}


//计算旋转矩阵2
Mat* Func_r(Mat* kk, Mat* p1,Mat* p2) {
    Mat k;
    if (kk->row != 1) {
        MatCopy(MatTrans(kk), &k);
    }
    else {
        MatCopy(kk, &k);
    }
    int count = 0;
    static Mat result;
    Mat  k_3, temp;
    MatCreate(&result, 1, 33);
    MatCreate(&k_3, 3, 3);
    MatCreate(&temp, 3, 3);
    MatSetVal(&k_3, k.element[0]);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++) {
            result.element[0][count] = k.element[0][3 * i + 0] * p1->element[0][j] + k.element[0][3 * i + 1] * p1->element[1][j] + k.element[0][3 * i + 2] * p1->element[2][j];
            count += 1;
        }
    MatCopy(MatMul(&k_3, MatTrans(&k_3)),&temp);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            result.element[0][count] = temp.element[i][j];
            count += 1;
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++) {
            result.element[0][count] = k.element[0][i + 0] * p2->element[0][j] + k.element[0][i + 3] * p2->element[1][j] + k.element[0][i + 6] * p2->element[2][j];
            count += 1;
        }
    MatDelete(&k);
    MatDelete(&k_3);
    MatDelete(&temp);
    return &result;
}
Mat* Deriv_r(Mat* abc,Mat* p1, Mat* p2, int n) {
    Mat x1, x2;
    Mat pp1, pp2;
    static Mat d;
    MatCopy(abc, &x1);
    MatCopy(abc, &x2);
    x1.element[n][0] -= 0.000001;
    x2.element[n][0] += 0.000001;
    MatCopy(Func_r(&x1, p1, p2),&pp1);
    MatCopy(Func_r(&x2, p1, p2),&pp2);
    MatCopy(MatTrans(MatNumMul(MatSub(&pp2, &pp1), 500000.0)),&d);
    MatDelete(&x1);
    MatDelete(&x2);
    MatDelete(&pp1);
    MatDelete(&pp2);
    return &d;
}
Mat* solve_r( Mat* p1,Mat* p2) {
    Mat y;
    MatCreate(&y, 1, 33);
    int count = 0;
    for (int i=0; i < 3; i++)
        for (int j=0; j < 4; j++) {
            y.element[0][count] = p2->element[i][j];
            count += 1;
        }
    for (int i = 0; i < 9; i++) {
        if (i % 4 == 0) {
            y.element[0][count] = 1;
        }
        else {
            y.element[0][count] = 0.0;
        }
        count += 1;
    }
    for (int i=0; i < 3; i++)
        for (int j=0; j < 4; j++) {
            y.element[0][count] = p1->element[i][j];
            count += 1;
        }

    // 设置迭代计数，总迭代数，方程个数
    int step = 0, conve = 20;
    double mse, mse_tmp, q;
    // 雅克比矩阵，误差矩阵
    Mat J, fx, fx_tmp, H, dx, xk_tmp;
    MatCreate(&J, 33, 9);
    MatCreate(&fx, 33, 9);
    MatCreate(&fx_tmp, 33, 1);
    MatCreate(&H, 3, 3);
    // 上一次迭代的均方差
    double lase_mse = 0.0;
    double u = 1, v = 2;
    //初始值设置
    Mat xk;
    Mat k;
    MatCreate(&xk, 9, 1);
    MatCreate(&k, 1, 9);
    double buff[] = { 1,0,0,0,1,0,0,0,1 };
    MatSetVal(&xk, buff);
    MatSetVal(&k, buff);
    while (conve) {
        mse = 0.0;
        mse_tmp = 0.0;
        step += 1;
        MatCopy(MatSub(Func_r(&k,p1,p2), &y),&fx);
        mse += MatSum(MatPow(&fx, 2));
        for (int j = 0; j < 9; j++) {
            MatColStead(&J, Deriv_r(&xk, p1,p2,j), j);
        }
        MatCopy(MatAdd(MatMul(MatTrans(&J), &J), MatNumMul(MatEye(9), u)),&H);
        MatCopy(MatMul(MatMul(MatInv(MatNumMul(&H, -1)), MatTrans(&J)), MatTrans(&fx)),&dx);
        MatCopy(&xk, &xk_tmp);
        MatCopy(MatAdd(&xk_tmp, &dx),&xk_tmp);
        MatCopy(MatSub(Func_r(&xk_tmp, p1, p2),&y),&fx_tmp);
        mse_tmp = MatSum(MatPow(MatCol(&fx_tmp, 0), 2));
        // 判断是否下降
        q = (mse - mse_tmp) / (MatNumMul(MatMul(MatTrans(&dx), MatSub(MatNumMul(&dx, u), MatMul(MatTrans(&J), MatTrans(&fx)))), 0.5)->element[0][0]);
        if (q > 0) {
            double s = 1.0 / 3;
            mse = mse_tmp;
            MatCopy(&xk_tmp, &xk);
            MatCopy(MatTrans(&xk), &k);
            double temp = 1.0 - (2 * q - 1) * (2 * q - 1) * (2 * q - 1);
            if (s > temp)
                u = u * s;
            else
                u = u * temp;
        }
        else {
            u = u * v;
            v = 2 * v;
            MatCopy(&xk, &xk_tmp);
            MatCopy(MatTrans(&xk), &k);
        }
        if (mse - lase_mse<0.000001 && mse - lase_mse>-0.000001) {
            static Mat result;
            MatCreate(&result, 3, 3);
            MatSetVal(&result, k.element[0]);

            MatDelete(&J);
            MatDelete(&fx);
            MatDelete(&fx_tmp);
            MatDelete(&H);
            MatDelete(&dx);
            MatDelete(&xk_tmp);
            MatDelete(&xk);
            MatDelete(&k);
            MatDelete(&y);
            return &result;
        }
        //printf("error:%f\n", (mse - lase_mse));
        lase_mse = mse;
        conve -= 1;
    }
    return NULL;
}


// 计算旋转角度
//构造函数
Mat* Func_rotation_angle(Mat* angle_kk) {
    Mat angle;
    if (angle_kk->row != 1) {
        MatCopy(MatTrans(angle_kk), &angle);
    }
    else {
        MatCopy(angle_kk, &angle);
    }
    double x, y, z;
    x = angle.element[0][0];
    y = angle.element[0][1];
    z = angle.element[0][2];
    static Mat result;
    MatCreate(&result, 1, 9);
    result.element[0][0] = cos(y) * cos(z) + sin(y) * sin(z);
    result.element[0][1] = -sin(z) * cos(y) + sin(x) * sin(y) * cos(z);
    result.element[0][2] = sin(y) * cos(x);

    result.element[0][3] = cos(x) * sin(z);
    result.element[0][4] = cos(x) * cos(z);
    result.element[0][5] = -sin(x);

    result.element[0][6] = -sin(y) * cos(z) + sin(x) * cos(y) * sin(z);
    result.element[0][7] = sin(y) * sin(z) + sin(x) * cos(y) * cos(z);
    result.element[0][8] = cos(x) * cos(y);
    MatDelete(&angle);
    return &result;
}
//计算微分
Mat* Deriv_rotation_angle(Mat* abc, int n) {
    Mat x1, x2;
    Mat p1, p2;
    static Mat d;
    MatCopy(abc, &x1);
    MatCopy(abc, &x2);
    x1.element[n][0] -= 0.000001;
    x2.element[n][0] += 0.000001;
    MatCopy(Func_rotation_angle(&x1),&p1);
    MatCopy(Func_rotation_angle(&x2),&p2);
    MatCopy(MatTrans(MatNumMul(MatSub(&p2, &p1), 500000.0)),&d);
    MatDelete(&x1);
    MatDelete(&x2);
    MatDelete(&p1);
    MatDelete(&p2);
    return &d;
}
//计算旋转角
Mat* solve_rotation_angle(Mat* pos, Mat* R, Mat* position, Mat* init) {
    // 去除位移后的姿态变化量
    Mat delt_att, r, new_p,trans_pos;
    MatCreate(&new_p, 4, 3);
    MatCreate(&delt_att, 4, 3);
    MatCreate(&r, 3, 3);
    MatCopy(delt_x(position), &delt_att);
    MatCopy(MatMul(MatInv(R), MatTrans(&delt_att)), &new_p);
    MatCopy(MatTrans(pos), &trans_pos);
    MatCopy(solve_r(&trans_pos, &new_p), &r);
    // 计算姿态角
    Mat y;
    MatCreate(&y, 1, 9);
    int count = 0;
    for (int i=0; i < 3; i++)
        for (int j=0; j < 3; j++) {
            y.element[0][count] = r.element[i][j];
            count += 1;
        }
    // 设置迭代计数，总迭代数，方程个数
    int step = 0, conve = 40;
    double mse, mse_tmp, q;
    // 雅克比矩阵，误差矩阵
    Mat J, fx, fx_tmp, H, dx, xk_tmp;
    MatCreate(&J, 9, 3);
    MatCreate(&fx, 9, 1);
    MatCreate(&fx_tmp, 9, 1);
    MatCreate(&H, 3, 3);
    // 上一次迭代的均方差
    double lase_mse = 0.0;
    double u = 1, v = 2;
    //初始值设置
    Mat xk;
    static Mat k;
    MatCreate(&xk, 3, 1);
    MatCreate(&k, 1, 3);
    double buff[] = { 0,0,0 };
    MatSetVal(&xk, buff);
    MatSetVal(&k, buff);
    while (conve) {
        mse = 0.0;
        mse_tmp = 0.0;
        step += 1;
        MatCopy(MatSub(Func_rotation_angle(&k), &y),&fx);
        mse += MatSum(MatPow(&fx, 2));
        for (int j = 0; j < 3; j++) {
            MatColStead(&J, Deriv_rotation_angle(&xk, j), j);
        }
        MatCopy(MatAdd(MatMul(MatTrans(&J), &J), MatNumMul(MatEye(3), u)),&H);
        MatCopy(MatMul(MatMul(MatInv(MatNumMul(&H, -1)), MatTrans(&J)), MatTrans(&fx)),&dx);
        MatCopy(&xk, &xk_tmp);
        MatCopy(MatAdd(&xk_tmp, &dx), &xk_tmp);
        MatCopy(MatSub(Func_rotation_angle(&xk_tmp), &y),&fx_tmp);
        mse_tmp = MatSum(MatPow(MatCol(&fx_tmp, 0), 2));
        // 判断是否下降
        q = (mse - mse_tmp) / (MatNumMul(MatMul(MatTrans(&dx), MatSub(MatNumMul(&dx, u), MatMul(MatTrans(&J), MatTrans(&fx)))), 0.5)->element[0][0]);
        if (q > 0) {
            double s = 1.0 / 3;
            mse = mse_tmp;
            MatCopy(&xk_tmp, &xk);
            MatCopy(MatTrans(&xk), &k);
            double temp = 1.0 - (2 * q - 1.0) * (2 * q - 1) * (2 * q - 1);
            if (s > temp)
                u = u * s;
            else
                u = u * temp;
        }
        else {
            u = u * v;
            v = 2 * v;
            MatCopy(&xk, &xk_tmp);
            MatCopy(MatTrans(&xk), &k);
        }
        if (mse - lase_mse<0.003 && mse - lase_mse>-0.003) {
            MatDelete(&delt_att);
            MatDelete(&r);
            MatDelete(&new_p);
            MatDelete(&trans_pos);
            MatDelete(&J);
            MatDelete(&fx);
            MatDelete(&fx_tmp);
            MatDelete(&H);
            MatDelete(&dx);
            MatDelete(&xk_tmp);
            MatDelete(&xk);
            MatDelete(&y);
            return &k;
        }
        printf("error:%f\n", (mse - lase_mse));
        lase_mse = mse;
        conve -= 1;   
    }
    return NULL;

}


// 异常点检测，暂时保留
int Outlier_detection(double* point){
    return 0;
}

// 计算当前位置
Mat* cal_position(Mat* angle,Mat* pos){
    int count = 0;
    Mat sol,y;
    static Mat position_list;
    MatCreate(&sol, 1, 4);
    MatCreate(&y, 1, 6);
    MatCopy(MatPow(sensor_distance(pos),2),&y);
    MatCreate(&position_list, 4, 3);
    // 基站到传感的单位向量
    Mat unit_ray, unit_ray_dot_list;
    MatCreate(&unit_ray, 4, 3);
    MatCreate(&unit_ray_dot_list, 1, 6);
    // 计算方向向量
    for (int i = 0; i < 4; i++)
        MatRowStead(&unit_ray, ray(angle->element[0][i], angle->element[0][4+i]), i);
    // 计算余弦值
    for(int i=0;i<4;i++)
        for (int j = i + 1; j < 4; j++) {
            unit_ray_dot_list.element[0][count] = MatDot(MatRow(&unit_ray, i), MatRow(&unit_ray, j));
            count += 1;
        }
    MatCopy(solve_x(&unit_ray_dot_list, &y),&sol);
    // 计算传感器坐标
    for (int i = 0; i < 4; i++) {
        MatRowStead(&position_list, MatNumMul(MatRow(&unit_ray, i), sol.element[0][i]), i);
    }
    MatDelete(&sol);
    MatDelete(&y);
    MatDelete(&unit_ray);
    MatDelete(&unit_ray_dot_list);
    return &position_list;

}

