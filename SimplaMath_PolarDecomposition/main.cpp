#include "SimpleMath.h"

#include <iostream>
#include <sstream>

void PrintMatrix(const char* _matName, const math::Matrix& _mat);
#define MATRIX_NAME(_matName) #_matName, _matName
//MatrixCmp(MATRIX_NAME(matA), MATRIX_NAME(matB));
void MatrixCmp(const char* _matNameA, const math::Matrix& _matA, const char* _matNameB, const math::Matrix& _matB, float _threshold = 0.01f);



int main()
{
	using namespace math;
	Matrix H
	{
		1.f, 0.f, 1.f, 0.f,
		1.f, 1.f, 1.f, 0.f,
		1.f, 1.f, 0.f, 1.f, 
		0.f, 0.f, 0.f, 1.f
	};
	H = H.Transpose();

	Matrix L = H;
	Matrix T = Matrix::CreateTranslation(H.Translation());

	L.Translation(Vector3::Zero);

	//H = LT
	MatrixCmp(MATRIX_NAME(H), MATRIX_NAME(L * T));
	
	Matrix R{};
	Matrix K{};

	L.DecomposePolar(R, K);
	//R = R.Transpose();
	//K = K.Transpose();

	if (R.Determinant() < 0.f)
	{
		R = -R;
		K = -K;
	}

	//L = KR
	//H = KRT
	MatrixCmp(MATRIX_NAME(L), MATRIX_NAME(K * R));
	MatrixCmp(MATRIX_NAME(H), MATRIX_NAME(K * R * T));


	return 0;
}

void PrintMatrix(const char* _matName, const math::Matrix& _mat)
{
	std::cout << "\n행렬 이름: " << _matName << std::endl;

	for (int i = 0; i < 4; ++i)
	{
		std::stringstream str{};

		for (int j = 0; j < 4; ++j)
		{
			str << std::to_string(_mat.m[i][j]) << ", ";
		}

		str << "\n";

		std::cout << str.str() << std::endl;
	}
	std::cout << "============================================\n\n\n" << std::endl;
}


void MatrixCmp(const char* _matNameA, const math::Matrix& _matA, const char* _matNameB, const math::Matrix& _matB, float _threshold)
{
	std::cout << "\n\n\n============================================" << std::endl;
	PrintMatrix(_matNameA, _matA);
	PrintMatrix(_matNameB, _matB);

	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			if (_threshold < fabsf(_matA.m[i][j] - _matB.m[i][j]))
			{
				std::cout << "=========================";
				std::cout << " DIFF ";
				std::cout << "=========================" << std::endl;
				return;
			}
		}
	}

	std::cout << "=========================";
	std::cout << " SAME ";
	std::cout << "=========================" << std::endl;
}

