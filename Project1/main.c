#include <stdio.h>
#include<math.h>
#include<time.h>
#include "light_matrix.h"
#include "location.h"
Mat pos, position, init_p[5];
// 测试使用Raw_data;
int* Raw_data;
// 主函数，接收传感器发送的时间信息，计算当前的位姿
void main(void) {
    // 测试专用
    ///F0 13 C2 F2 13 AF F4 14 9E F6 14 99 F1 10 67 F3 0F F1 F5 0F F0 F7 10 68
    int Raw_data[96] = { 0xf0,0x13,0xc2,0xf2,0x13,0xaf,0xf4,0x14,0x9e,0xf6,0x14,0x99,0xf1,0x10,0x67,0xf3,0x0F,0xf1,0xf5,0x0F,0xF0,0xf7,0x10,0x68,
                         0xf0,0x13,0xc2,0xf2,0x13,0xaf,0xf4,0x14,0x9e,0xf6,0x14,0x99,0xf1,0x10,0x67,0xf3,0x0F,0xf1,0xf5,0x0F,0xF0,0xf7,0x10,0x68,
                         0xf0,0x13,0xc2,0xf2,0x13,0xaf,0xf4,0x14,0x9e,0xf6,0x14,0x99,0xf1,0x10,0x67,0xf3,0x0F,0xf1,0xf5,0x0F,0xF0,0xf7,0x10,0x68,
                         0xf0,0x13,0xc2,0xf2,0x13,0xaf,0xf4,0x14,0x9e,0xf6,0x14,0x99,0xf1,0x10,0x67,0xf3,0x0F,0xf1,0xf5,0x0F,0xF0,0xf7,0x10,0x68,
    };
    // 计时测试
    int start, end;
    int nums = sizeof(Raw_data) / 4;
    Mat time_mean, init, R, rotation_angle;
    MatCreate(&R, 3, 3);
    MatCreate(&init, 4, 3);
    MatCreate(&rotation_angle, 1, 3);
    int count = -1;
    MatCreate(&time_mean, 1, 8);
    MatCreate(&pos, 4, 3);
    MatCreate(&position, 4, 1);
    MatCreate(&init, 4, 3);
    double pol_val[] = { -36.5,-36.5,0,
                         -36.5, 36.5,0,
                          36.5, 36.5,0,
                          36.5,-36.5,0, };
    MatSetVal(&pos, pol_val);
    while (1) {

        start = clock();
        // 数据预处理
        Mat data = *filter(Raw_data, nums);
        // 计算角度信息
        tick2angle(&data);
        // 计算角度的平均值
        time_mean = *cal_mean(&data);
        // 计算位置
        position = *cal_position(&time_mean,&pos);
        count++;
        // 计算转换矩阵
        if (count < 5) {
            MatCopy(&position, &init_p[count]);
        }
        else if (count < 6) {
            MatCopy(cal_init_mean(init_p), &init);
            // 计算旋转矩阵
            MatCopy(solve_R(&position, &pos), &R);
        }
        else {
            printf("--------------\n");
            MatCopy(solve_rotation_angle(&pos, &R, &position, &init), &rotation_angle);
            printf("Position_Mat\n");
            MatDump(&position);
            printf("-----------------------------------------------------\n");
            printf("Rotation_Angle_Mat\n");
            MatDump(&rotation_angle);
            printf("-----------------------------------------------------\n");
            end = clock();
            printf("time=%d ms\n", end - start);
            printf("-----------------------------------------------------\n");
        }
    }
    
}