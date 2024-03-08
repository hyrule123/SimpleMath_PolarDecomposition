#include "SimpleMath.h"

#include <iostream>
#include <sstream>

#include "Eigen/Dense"

void PrintMatrix(const char* _matName, const math::Matrix& _mat);

#define MATRIX_NAME(_matName) #_matName, _matName
//MatrixCmp(MATRIX_NAME(matA), MATRIX_NAME(matB));
void MatrixCmp(const char* _matNameA, const math::Matrix& _matA, const char* _matNameB, const math::Matrix& _matB, float _threshold = 0.01f);


//https://math.stackexchange.com/questions/237369/given-this-transformation-matrix-how-do-i-decompose-it-into-translation-rotati/3554913#3554913
int main()
{
	//rotation matrix Original
	//float sqrt2Inv = 1.f / sqrt2;
	//math::Matrix answerRotationMatrix =
	//{
	//	sqrt2Inv, sqrt2Inv, 0.f, 0.f,
	//	-sqrt2Inv, sqrt2Inv, 0.f, 0.f,
	//	0.f, 0.f, 1.f, 0.f,
	//	0.f, 0.f, 0.f, 1.f
	//};
	////->{ 0.f, 0.f, 0.785398245f }


	float sqrt2 = ::sqrtf(2.f);
	math::Vector3 answerScale{ -2.f * sqrt2, sqrt2, 1.f };
	math::Vector3 answerRot{ 0.f, 0.f, 0.785398245f };
	math::Vector3 answerTranslation{ 2.f, -1.f, 0.f };

	math::Matrix answer =
		math::Matrix::CreateScale(answerScale)
		* math::Matrix::CreateFromYawPitchRoll(answerRot)
		* math::Matrix::CreateTranslation(answerTranslation);

	//Sign test
	math::Matrix test = answer * answer.Transpose();
	float scaleSign[3];
	for (int i = 0; i < 3; ++i)
	{
		scaleSign[i] = ::sqrtf(test.m[i][i]);
		if (std::signbit(answer.m[i][i]))
		{
			scaleSign[i] *= -1.f;
		}
	}

	//Matrix of multiple transformation matrix being multiplied
	//math::Matrix childWorld{};
	//{
	//	using namespace math;
	//	Matrix child = Matrix::Identity;
	//	Matrix parent = Matrix::Identity;

	//	Vector3 parentScale = Vector3(1.f, 2.f, 3.f);
	//	Quaternion parentRot = Quaternion::CreateFromYawPitchRoll(Vector3(XM_PI / 1.f, XM_PI / 2.f, XM_PI / 3.f));
	//	Vector3 parentPos = Vector3(3.f, 2.f, 1.f);

	//	parent *= Matrix::CreateScale(parentScale);
	//	parent *= Matrix::CreateFromQuaternion(parentRot);
	//	parent *= Matrix::CreateTranslation(parentPos);

	//	Vector3 childScale = Vector3(4.f, 3.f, 2.f);
	//	Quaternion childRot = Quaternion::CreateFromYawPitchRoll(Vector3(XM_PI * 3.1235f, XM_PI * 12.3853f, XM_PI * 8.192f));
	//	Vector3 childPos = Vector3(4.f, 8.f, 12.f);
	//	child *= Matrix::CreateScale(childScale);
	//	child *= Matrix::CreateFromQuaternion(childRot);
	//	child *= Matrix::CreateTranslation(childPos);

	//	childWorld = child * parent;
	//}


	//================================ Decompose procedure =================================================
	//Set 'H' matrix what you want to decompose.
	math::Matrix H = answer;
	
	math::Matrix L = H;
	math::Matrix T = math::Matrix::CreateTranslation(H.Translation());

	//Position
	math::Vector3 pos = T.Translation();

	L.Translation(math::Vector3::Zero);

	//H = LT
	MatrixCmp(MATRIX_NAME(H), MATRIX_NAME(L * T));
	
	math::Matrix R{};
	math::Matrix K{};

	L.DecomposePolar(R, K);


	if (R.Determinant() < 0.f)
	{
		for (int r = 0; r < 3; ++r)
		{
			for (int c = 0; c < 3; ++c)
			{
				R.m[r][c] = -(R.m[r][c]);
				K.m[r][c] = -(K.m[r][c]);
			}
		}
	}


	//L = KR
	//H = KRT
	MatrixCmp(MATRIX_NAME(L), MATRIX_NAME(K * R));
	MatrixCmp(MATRIX_NAME(H), MATRIX_NAME(K * R * T));

	Eigen::Matrix<float, 4, 4, Eigen::RowMajor> eigMat{};
	eigMat = Eigen::Map<Eigen::Matrix<float, 4, 4, Eigen::RowMajor>>(K.m[0], 4, 4);
	
	Eigen::EigenSolver<Eigen::Matrix4Xf> solver(eigMat);

	assert(solver.info() == Eigen::Success);

	auto eigVec = solver.eigenvectors();
	auto eigVal = solver.eigenvalues();

	//may be scale value?
	math::Vector4 f{};
	f.x = eigVal(0).real();
	f.y = eigVal(1).real();
	f.z = eigVal(2).real();
	f.w = eigVal(3).real();
	math::Vector3 scale = { f.x, f.y, f.z };

	math::Matrix X{};
	for (int r = 0; r < 4; ++r)
	{
		for (int c = 0; c < 4; ++c)
		{
			X.m[r][c] = eigVec(r * 4 + c).real();
		}
	}

	
	for (int i = 0; i < 3; ++i)
	{
		float* pAnswerScaleAxis{};
		float* pScaleAxis{};
		switch (i)
		{
		case 0:
			pScaleAxis = &(scale.x);
			pAnswerScaleAxis = &(answerScale.x);
			break;
		case 1:
			pScaleAxis = &(scale.y);
			pAnswerScaleAxis = &(answerScale.y);
			break;
		case 2:
			pScaleAxis = &(scale.z);
			pAnswerScaleAxis = &(answerScale.z);
			break;
		default:
			assert(false);
			break;
		}

		float& answerScaleAxis = *pAnswerScaleAxis;
		float& scaleAxis = *pScaleAxis;

		if (0.f > answerScaleAxis * scaleAxis)
		{
			scaleAxis *= -1.f;
			for (int j = 0; j < 3; ++j)
			{
				R.m[i][j] *= -1.f;
			}
		}
	}


	//Combined transformation matrix cannot use this decomposition method.
	math::Matrix recomposed =
		math::Matrix::CreateScale(scale)
		* R
		* math::Matrix::CreateTranslation(pos);

	MatrixCmp(MATRIX_NAME(answer), MATRIX_NAME(recomposed));



	return 0;
}

void PrintMatrix(const char* _matName, const math::Matrix& _mat)
{
	std::cout << "\nMATRIX NAME: " << _matName << std::endl;

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

