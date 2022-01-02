#include "light_matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define MAT_LEGAL_CHECKING

//#define min(a, b) ((a) > (b) ? (b) : (a))
#define equal(a, b)	((a-b)<1e-7 && (a-b)>-(1e-7))

/************************************************************************/
/*                          Private Function                            */
/************************************************************************/

void swap(int *a, int *b)
{
	int m;
	m = *a;
	*a = *b;
	*b = m;
}
 
void perm(int list[], int k, int m, int* p, Mat* mat, double* det) 
{
	int i;

	if(k > m){
		double res = mat->element[0][list[0]];

		for(i = 1; i < mat->row ; i++){
			res *= mat->element[i][list[i]];
		}

		if(*p%2){
			//odd is negative
			*det -= res;
		}else{
			//even is positive
			*det += res;
		}
	}
	else{
		// if the element is 0, we don't need to calculate the value for this permutation
		if(!equal(mat->element[k][list[k]], 0.0f))
			perm(list, k + 1, m, p, mat, det);
		for(i = k+1; i <= m; i++)
		{
			if(equal(mat->element[k][list[i]], 0.0f))
				continue;
			swap(&list[k], &list[i]);
			*p += 1;
			perm(list, k + 1, m, p, mat, det);
			swap(&list[k], &list[i]);
			*p -= 1; 
		}
	}
}

/************************************************************************/
/*                           Public Function                            */
/************************************************************************/

Mat* MatCreate(Mat* mat, int row, int col)
{
	if (mat->flag == 1) {
		MatDelete(mat);
	}
	int i;
	mat->element = (double**)malloc(row * sizeof(double*));
	if(mat->element == NULL){
		printf("mat create fail!\n");
		return NULL;
	}
	for(i = 0 ; i < row ; i++){
		mat->element[i] = (double*)malloc(col * sizeof(double));	
		if(mat->element[i] == NULL){
			int j;
			printf("mat create fail!\n");
			for(j = 0 ; j < i ; j++)
				free(mat->element[j]);
			free(mat->element);
			return NULL;
		}
	}

	mat->row = row;
	mat->col = col;
	mat->flag = 1;

	return mat;
}

void MatDelete(Mat* mat)
{
	int i;
	mat->flag = 0;
	for(i = 0 ; i<mat->row ; i++)
		free(mat->element[i]);
	free(mat->element);
}

Mat* MatSetVal(Mat* mat, double* val)
{
	int row,col;

	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
			mat->element[row][col] = val[col + row * mat->col];
		}
	}

	return mat;
}

void MatDump(const Mat* mat)
{
	int row,col;

#ifdef MAT_LEGAL_CHECKING
	if(mat == NULL){
		return ;
	}
#endif

	printf("Mat %dx%d:\n", mat->row, mat->col);
	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
			printf("%.12f\t", mat->element[row][col]);
		}
		printf("\n");
	}
}

Mat* MatZeros(Mat* mat)
{
	int row,col;
	for(row = 0 ; row < mat->row ; row++){
		for(col = 0 ; col < mat->col ; col++){
			mat->element[row][col] = 0.0f;
		}
	}

	return mat;
}

Mat* MatEye(int k)
{
	int i;
	static Mat mat;
	MatCreate(&mat, k, k);
	MatZeros(&mat);
	for(i = 0 ; i < k ; i++){
		mat.element[i][i] = 1.0f;
	}

	return &mat;
}

/* dst = src1 + src2 */
static int add_flag = 0;
Mat* MatAdd(Mat* src1, Mat* src2)
{
	int row, col;
	Mat dst;
	static Mat dst0,dst1;
	MatCreate(&dst, src1->row, src1->col);

#ifdef MAT_LEGAL_CHECKING
	if( !(src1->row == src2->row && src2->row == dst.row && src1->col == src2->col && src2->col == dst.col) ){
		printf("err check, unmatch matrix for MatAdd\n");
		MatDump(src1);
		MatDump(src2);
		MatDump(&dst);
		MatDelete(&dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < src1->row ; row++){
		for(col = 0 ; col < src1->col ; col++){
			dst.element[row][col] = src1->element[row][col] + src2->element[row][col];
		}
	}
	if (add_flag== 0) {
		MatCopy(&dst, &dst0);
		add_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		add_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

/* dst = src1 - src2 */
static int sub_flag = 0;
Mat* MatSub(Mat* src1, Mat* src2)
{
	int row, col;
	Mat dst;
	static Mat dst0, dst1;
	MatCreate(&dst, src1->row, src1->col);

#ifdef MAT_LEGAL_CHECKING
	if( !(src1->row == src2->row && src2->row == dst.row && src1->col == src2->col && src2->col == dst.col) ){
		printf("err check, unmatch matrix for MatSub\n");
		MatDump(src1);
		MatDump(src2);
		MatDump(&dst);
		MatDelete(&dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < src1->row ; row++){
		for(col = 0 ; col < src1->col ; col++){
			dst.element[row][col] = src1->element[row][col] - src2->element[row][col];
		}
	}
	if (sub_flag == 0) {
		MatCopy(&dst, &dst0);
		sub_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		sub_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
	
}

/* dst = src1 * src2 */
static int mul_flag = 0;
Mat* MatMul(Mat* src1, Mat* src2)
{
	int row, col;
	int i;
	double temp;
	Mat dst;
	static Mat dst0, dst1;
	MatCreate(&dst, src1->row, src2->col);

#ifdef MAT_LEGAL_CHECKING
	if( src1->col != src2->row || src1->row != dst.row || src2->col != dst.col ){
		printf("err check, unmatch matrix for MatMul\n");
		MatDump(src1);
		MatDump(src2);
		MatDump(&dst);
		MatDelete(&dst);
		return NULL;
	}
#endif

	for(row = 0 ; row < dst.row ; row++){
		for(col = 0 ; col < dst.col ; col++){
			temp = 0.0f;
			for(i = 0 ; i < src1->col ; i++){
				temp += src1->element[row][i] * src2->element[i][col];
			}
			dst.element[row][col] = temp;
		}
	}
	if (mul_flag == 0) {
		MatCopy(&dst, &dst0);
		mul_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		mul_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

/* dst = src' */
static int trans_flag = 0;
Mat* MatTrans(Mat* src)
{
	int row, col;
	Mat dst;
	static Mat dst0, dst1;
	MatCreate(&dst, src->col, src->row);
	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++){
			dst.element[col][row] = src->element[row][col];
		}
	}
	if (trans_flag == 0) {
		MatCopy(&dst, &dst0);
		trans_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		trans_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

// return det(mat)
double MatDet(Mat* mat)
{
	double det = 0.0f;
	int plarity = 0;
	int *list;
	int i;

#ifdef MAT_LEGAL_CHECKING
	if( mat->row != mat->col){
		printf("err check, not a square matrix for MatDetermine\n");
		MatDump(mat);
		return 0.0f;
	}
#endif

	list = (int*)malloc(sizeof(int)*mat->col);
	if(list == NULL){
		printf("malloc list fail\n");
		return -1;
	}
	for(i = 0 ; i < mat->col ; i++)
		list[i] = i;

	perm(list, 0, mat->row-1, &plarity, mat, &det);
	free(list);

	return det;
}

// dst = adj(src)
Mat* MatAdj(Mat* src)
{
	Mat smat;
	int row, col;
	int i,j,r,c;
	double det;
	static Mat dst;
	MatCreate(&dst, src->row, src->col);

#ifdef MAT_LEGAL_CHECKING
	if( src->row != src->col || src->row != dst.row || src->col != dst.col){
		printf("err check, not a square matrix for MatAdj\n");
		MatDump(src);
		MatDump(&dst);
		return NULL;
	}
#endif

	MatCreate(&smat, src->row-1, src->col-1);

	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++){
			r = 0;
			for(i = 0 ; i < src->row ; i++){
				if(i == row)
					continue;
				c = 0;
				for(j = 0; j < src->col ; j++){
					if(j == col)
						continue;
					smat.element[r][c] = src->element[i][j];
					c++;
				}
				r++;
			}
			det = MatDet(&smat);
			if((row+col)%2)
				det = -det;
			dst.element[col][row] = det;
		}
	}

	MatDelete(&smat);

	return &dst;
}

// dst = src^(-1)
Mat* MatInv(Mat* src)
{
	Mat adj_mat;
	double det;
	int row, col;
	static Mat dst;
	MatCreate(&dst, src->row, src->col);

#ifdef MAT_LEGAL_CHECKING
	if( src->row != src->col || src->row != dst.row || src->col != dst.col){
		printf("err check, not a square matrix for MatInv\n");
		MatDump(src);
		MatDump(&dst);
		return NULL;
	}
#endif
	MatCreate(&adj_mat, src->row, src->col);
	MatCopy(MatAdj(src), &adj_mat);
	det = MatDet(src);

	if(equal(det, 0.0f)){
		printf("err, determinate is 0 for MatInv\n");
		return NULL;
	}
	
	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++)
			dst.element[row][col] = adj_mat.element[row][col]/det;
	}
	
	MatDelete(&adj_mat);

	return &dst;
}

void MatCopy(Mat* src, Mat* dst)
{
	int row, col;
	if (dst->flag == 1) {
		MatDelete(dst);
	}
	
	MatCreate(dst, src->row, src->col);
	for(row = 0 ; row < src->row ; row++){
		for(col = 0 ; col < src->col ; col++)
			dst->element[row][col] = src->element[row][col];
	}
}

//之后添加的一些函数
// 数乘运算
Mat* MatNumMul(Mat* src1, double k)
{
	int row, col;
	static Mat dst;
	MatCreate(&dst, src1->row, src1->col);
	for (row = 0; row < dst.row; row++) {
		for (col = 0; col < dst.col; col++) {
			dst.element[row][col] = src1->element[row][col] * k;
		}
	}
	return &dst;
}

// 向量单位化
Mat* MatUnit(Mat* src1)
{
	int row, col;
	static Mat dst;
	MatCreate(&dst, src1->row, src1->col);

#ifdef MAT_LEGAL_CHECKING
	if (src1->col != 1 && src1->row != 1) {
		printf("err check, unmatch matrix for MatUnit\n");
		MatDump(src1);
		MatDump(&dst);
		return NULL;
	}
#endif
	double delt;
	delt = sqrt(MatSum(MatPow(src1, 2))) ;
	for (row = 0; row < dst.row; row++) {
		for (col = 0; col < dst.col; col++) {
			dst.element[row][col] = src1->element[row][col] / delt;
		}
	}
	return &dst;
}

// 矩阵行堆叠
Mat* MatRowStock(Mat* src1, Mat* src2)
{
	int row, col;
	static Mat dst;
	MatCreate(&dst, src1->row+src2->row,src1->col);

#ifdef MAT_LEGAL_CHECKING
	if (src1->col != src2->col) {
		printf("err check, unmatch matrix for MatRowStock\n");
		MatDump(src1);
		MatDump(src2);
		return NULL;
	}
#endif
	for(row=0;row<src1->row;row++)
		for (col = 0; col < src1->col; col++) {
			dst.element[row][col] = src1->element[row][col] ;
		}
	for (row = src1->row; row < (src1->row+ src2->row); row++)
		for (col = 0; col < src1->col; col++) {
			dst.element[row][col] = src2->element[row-src1->row][col];
		}

	return &dst;
}

// 矩阵列堆叠
Mat* MatColStock(Mat* src1, Mat* src2)
{
	int row, col;
	static Mat dst;
	MatCreate(&dst, src1->row , src1->col+src2->col);

#ifdef MAT_LEGAL_CHECKING
	if (src1->row != src2->row) {
		printf("err check, unmatch matrix for MatColStock\n");
		MatDump(src1);
		MatDump(src2);
		return NULL;
	}
#endif
	for (col = 0; col < src1->col; col++)
		for (row = 0; row < src1->row; row++) {
			dst.element[row][col] = src1->element[row][col];
		}
	for (col = src1->col; col < (src1->col + src2->col); col++)
		for (row = 0; row < src1->row; row++) {
			dst.element[row][col] = src2->element[row][col-src1->col];
		}
	return &dst;
}

// 所有元素求和
double MatSum(Mat* src1)
{
	int row, col;
	double result = 0.0;
	for (row = 0; row < src1->row; row++) {
		for (col = 0; col < src1->col; col++) {
			result += src1->element[row][col] ;
		}
	}
	return result;
}

// 次方运算
static int pow_flag = 0;
Mat* MatPow(Mat* src1, double k)
{
	int row, col;
	Mat dst;
	static Mat dst0, dst1;
	MatCreate(&dst, src1->row, src1->col);
	for (row = 0; row < src1->row; row++) {
		for (col = 0; col < src1->col; col++) {
			dst.element[row][col] = pow(src1->element[row][col],k) ;
		}
	}
	if (pow_flag == 0) {
		MatCopy(&dst, &dst0);
		pow_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		pow_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

// 矩阵行替换
void MatRowStead(Mat* src1, Mat* src2,int k)
{
	int col;

#ifdef MAT_LEGAL_CHECKING
	if (src1->col != src2->col||src2->row != 1||k>=src1->row) {
		printf("err check, unmatch matrix for MatRowStead\n");
		MatDump(src1);
		MatDump(src2);
	}
#endif
	for (col = 0; col < src1->col; col++) {
			src1->element[k][col] = src2->element[0][col];
	}

}

// 矩阵列替换
void MatColStead(Mat* src1, Mat* src2, int k)
{
	int row;


#ifdef MAT_LEGAL_CHECKING
	if (src1->row != src2->row || src2->col != 1 || k >= src1->col) {
		printf("err check, unmatch matrix for MatColStead\n");
		MatDump(src1);
		MatDump(src2);
	}
#endif
	for (row = 0; row < src1->row; row++) {
		src1->element[row][k] = src2->element[row][0];
	}

}

//两个矩阵数乘
static mul2_flag = 0;
Mat* MatMul2(Mat* src1, Mat* src2)
{
	int row, col;
	Mat dst;
	static Mat dst0,dst1;
	MatCreate(&dst, src1->row, src1->col);

#ifdef MAT_LEGAL_CHECKING
	if (!(src1->row == src2->row && src2->row == dst.row && src1->col == src2->col && src2->col == dst.col)) {
		printf("err check, unmatch matrix for MatMul2\n");
		MatDump(src1);
		MatDump(src2);
		return NULL;
	}
#endif

	for (row = 0; row < src1->row; row++) {
		for (col = 0; col < src1->col; col++) {
			dst.element[row][col] = src1->element[row][col] * src2->element[row][col];
		}
	}
	if (mul2_flag == 0) {
		MatCopy(&dst, &dst0);
		mul2_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		mul2_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

// 矩阵取单独一行
static int row_flag = 0;
Mat* MatRow(Mat* src1, int k)
{
	int col;
	Mat dst;
	static Mat dst0,dst1;
	MatCreate(&dst, 1, src1->col);

#ifdef MAT_LEGAL_CHECKING
	if (k >= src1->row) {
		printf("err check, unmatch matrix for MatRow\n");
		MatDump(src1);
		return NULL;
	}
#endif
	for (col = 0; col < src1->col; col++) {
		dst.element[0][col] = src1->element[k][col];
	}
	if (row_flag == 0) {
		MatCopy(&dst, &dst0);
		row_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		row_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

// 矩阵取单独一列
static col_flag = 0;
Mat* MatCol(Mat* src1, int k)
{
	int row;
	Mat dst;
	static Mat dst0, dst1;
	MatCreate(&dst, src1->row, 1);

#ifdef MAT_LEGAL_CHECKING
	if (k >= src1->col) {
		printf("err check, unmatch matrix for MatCol\n");
		MatDump(src1);
		MatDelete(&dst);
		return NULL;
	}
#endif
	for (row = 0; row < src1->row; row++) {
		dst.element[row][0] = src1->element[row][k];
	}
	if (col_flag == 0) {
		MatCopy(&dst, &dst0);
		col_flag = 1;
		MatDelete(&dst);
		return &dst0;
	}
	else {
		MatCopy(&dst, &dst1);
		col_flag = 0;
		MatDelete(&dst);
		return &dst1;
	}
}

// 3维向量叉乘,默认是横着的
Mat* MatCross3(Mat* src1, Mat* src2)
{
	static Mat dst;
	MatCreate(&dst, 1,3);

#ifdef MAT_LEGAL_CHECKING
	if (src1->col != 3 || src2->col!=3 || src1->row!=1||src2->row != 1) {
		printf("err check, unmatch matrix for MatCross3\n");
		MatDump(src1);
		return NULL;
	}
#endif
	dst.element[0][0] = src1->element[0][1] * src2->element[0][2]- src1->element[0][2] * src2->element[0][1];
	dst.element[0][1] = -src1->element[0][0] * src2->element[0][2]+ src1->element[0][2] * src2->element[0][0];
	dst.element[0][2] = src1->element[0][0] * src2->element[0][1]- src1->element[0][1] * src2->element[0][0];

	return &dst;
}

// 向量点乘
double MatDot(Mat* src1, Mat* src2) {
	double dst=0.0;
#ifdef MAT_LEGAL_CHECKING
	if (src1->col != src2->col || src1->row != 1 || src2->row != 1) {
		printf("err check, unmatch matrix for MatDot\n");
		MatDump(src1);
		MatDump(src2);
		return -1;
	}
#endif
	for (int i = 0; i < src1->col; i++)
		dst += src1->element[0][i] * src2->element[0][i];
		return dst;
}

//计算每一列的平均值
Mat* MatColMean(Mat* src1) {
	static Mat dst;
	MatCreate(&dst, 1, src1->col);
	for (int i = 0; i < src1->col; i++) {
		dst.element[0][i] = MatSum(MatCol(src1, i)) / src1->row;
	}
	return &dst;
}