#ifndef __LIGHT_MATRIX__
#define __LIGHT_MATRIX__

typedef struct  {
	int row, col;
	int flag;
	double **element;
}Mat;

Mat* MatCreate(Mat* mat, int row, int col);
void MatDelete(Mat* mat);
Mat* MatSetVal(Mat* mat, double* val);
void MatDump(const Mat* mat);

Mat* MatZeros(Mat* src);
Mat* MatEye(int k);

Mat* MatAdd(Mat* src1, Mat* src2);
Mat* MatSub(Mat* src1, Mat* src2);
Mat* MatMul(Mat* src1, Mat* src2);
Mat* MatMul2(Mat* src1, Mat* src2);
Mat* MatNumMul(Mat* src1, double k);
Mat* MatPow(Mat* src1, double k);
Mat* MatUnit(Mat* src1);
Mat* MatTrans(Mat* src);
Mat* MatColStock(Mat* src1, Mat* src2);
Mat* MatRowStock(Mat* src1, Mat* src2);
void MatRowStead(Mat* src1, Mat* src2, int k);
void MatColStead(Mat* src1, Mat* src2, int k);
Mat* MatRow(Mat* src1, int k);
Mat* MatCol(Mat* src1, int k);
Mat* MatCross3(Mat* src1, Mat* src2);
Mat* MatColMean(Mat* src);
double MatDot(Mat* src1, Mat* src2);
double MatDet(Mat* mat);
double MatSum(Mat* src1);
Mat* MatAdj(Mat* src);
Mat* MatInv(Mat* src);

void MatCopy(Mat* src, Mat* dst);

#endif
