#ifndef _PLANE_FITTING_H
#define _PLANE_FITTING_H

/*
bool svd(vector<vector<float> > A, int K, vector<vector<float> > &U, vector<float> &S, vector<vector<float> > &V);

A: 输入待分解矩阵
K: 输入，取前K大奇异值及奇异向量
U[0],U[1],...,U[K-1]: 前K大奇异值对应的左奇异向量
S[0],S[1],...,S[K-1]: 前K大奇异值 S[0]>=S[1]>=...>=S[K-1]
V[0],V[1],...,V[K-1]: 前K大奇异值对应的右奇异向量
*/
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>

#include "tpoint.h"

namespace CPU_SVD {
	using namespace std;
	const int c_max_iter = 100000;
	const float c_eps = 0.0000001;

	float get_norm(float *x, int n) {
		float r = 0;
		for (int i = 0; i < n; i++)
			r += x[i] * x[i];
		return sqrt(r);
	}

	float normalize(float *x, int n) {
		float r = get_norm(x, n);
		if (r < c_eps)
			return 0;
		for (int i = 0; i < n; i++)
			x[i] /= r;
		return r;
	}

	inline float product(float*a, float *b, int n) {
		float r = 0;
		for (int i = 0; i < n; i++)
			r += a[i] * b[i];
		return r;
	}

	void orth(float *a, float *b, int n) {//|a|=1
		float r = product(a, b, n);
		for (int i = 0; i < n; i++)
			b[i] -= r*a[i];

	}

	template<typename _Tp>
	int calcCovarMatrix(const vector<vector<_Tp>>& mat, vector<vector<_Tp>>& covar, vector<_Tp>& mean, bool scale = false)
	{
		const int rows = mat.size();
		const int cols = mat[0].size();
		const int nsamples = rows;
		double scale_ = 1.;
		if (scale) scale_ = 1. / (nsamples /*- 1*/);

		covar.resize(cols);
		for (int i = 0; i < cols; ++i)
			covar[i].resize(cols, (_Tp)0);
		mean.resize(cols, (_Tp)0);

		for (int w = 0; w < cols; ++w) {
			for (int h = 0; h < rows; ++h) {
				mean[w] += mat[h][w];
			}
		}

		for (auto& value : mean) {
			value = 1. / rows * value;
		}

		for (int i = 0; i < cols; ++i) {
			vector<_Tp> col_buf(rows, (_Tp)0);
			for (int k = 0; k < rows; ++k)
				col_buf[k] = mat[k][i] - mean[i];

			for (int j = 0; j < cols; ++j) {
				double s0 = 0;
				for (int k = 0; k < rows; ++k) {
					s0 += col_buf[k] * (mat[k][j] - mean[j]);
				}
				covar[i][j] = (_Tp)(s0 * scale_);
			}
		}

		return 0;
	}


	bool svd(vector<vector<float> > A, int K, vector<vector<float> > &U, vector<float> &S, vector<vector<float> > &V) {
		int M = A.size();
		int N = A[0].size();
		U.clear();
		V.clear();
		S.clear();
		S.resize(K, 0);
		U.resize(K);
		for (int i = 0; i < K; i++)
			U[i].resize(M, 0);
		V.resize(K);
		for (int i = 0; i < K; i++)
			V[i].resize(N, 0);


		srand(time(0));
		float *left_vector = new float[M];
		float *next_left_vector = new float[M];
		float *right_vector = new float[N];
		float *next_right_vector = new float[N];
		int col = 0;
		for (int col = 0; col < K; col++) {
			float diff = 1;
			float r = -1;
			while (1) {
				for (int i = 0; i < M; i++)
					left_vector[i] = (float)rand() / RAND_MAX;
				if (normalize(left_vector, M) > c_eps)
					break;
			}

			for (int iter = 0; diff >= c_eps && iter < c_max_iter; iter++) {
				memset(next_left_vector, 0, sizeof(float)*M);
				memset(next_right_vector, 0, sizeof(float)*N);
				for (int i = 0; i < M; i++)
					for (int j = 0; j < N; j++)
						next_right_vector[j] += left_vector[i] * A[i][j];

				r = normalize(next_right_vector, N);
				if (r < c_eps) break;
				for (int i = 0; i < col; i++)
					orth(&V[i][0], next_right_vector, N);
				normalize(next_right_vector, N);

				for (int i = 0; i < M; i++)
					for (int j = 0; j < N; j++)
						next_left_vector[i] += next_right_vector[j] * A[i][j];
				r = normalize(next_left_vector, M);
				if (r < c_eps) break;
				for (int i = 0; i < col; i++)
					orth(&U[i][0], next_left_vector, M);
				normalize(next_left_vector, M);
				diff = 0;
				for (int i = 0; i < M; i++) {
					float d = next_left_vector[i] - left_vector[i];
					diff += d*d;
				}

				memcpy(left_vector, next_left_vector, sizeof(float)*M);
				memcpy(right_vector, next_right_vector, sizeof(float)*N);
			}
			if (r >= c_eps) {
				S[col] = r;
				memcpy((char *)&U[col][0], left_vector, sizeof(float)*M);
				memcpy((char *)&V[col][0], right_vector, sizeof(float)*N);
			}
			else {
				break;
			}
		}
		delete[] next_left_vector;
		delete[] next_right_vector;
		delete[] left_vector;
		delete[] right_vector;

		return true;
	}

} // namespace CPU_SVD

void PlaneFitting(std::vector<std::vector<float>> pt_matrix, Point3f &plane_norm)  //输入点的矩阵 按行存储
{
	std::vector<std::vector<float>> covar;   //协方差阵
	std::vector<float> mean;   //均值

	CPU_SVD::calcCovarMatrix(pt_matrix, covar, mean);

	std::vector<std::vector<float> > U;
	std::vector<float> S;
	std::vector<std::vector<float> > V;

	CPU_SVD::svd(covar,3,U,S,V);   //因为输入是三维点 协方差阵为 3*3 的方阵

	plane_norm.x = U[2][0];
	plane_norm.y = U[2][1];
	plane_norm.z = U[2][2];//输出的奇异值 按照从大到小排列 直接取最后一个

}

#endif // _PLANE_FITTING_H