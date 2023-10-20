#include<iostream>
#include<vector>

using namespace std;

class Matrix_and_fun {
private:
	vector<vector<double>> a;
	vector<double> b;
	vector<vector<double>> x_GE;
	vector<double> x_LU;
	vector<double> y_LU;
	vector<vector<double>> L;
	vector<vector<double>> U;
	vector<double> GE;
public:
	void setA(vector<vector<double>> matrix) {
		this->a = matrix;
	}

	vector<vector<double>> getA() {
		return this->a;
	}

	void setB(vector<double> matrix) {
		this->b = matrix;
	}

	vector<double> getB() {
		return this->b;
	}

	vector<double> getGE() {
		return this->GE;
	}

	vector<double> get_yLU() {
		return this->y_LU;
	}

	vector<double> get_xLU() {
		return this->x_LU;
	}

	vector<vector<double>> getX_ge() {
		return this->x_GE;
	}

	vector<vector<double>> getL() {
		return this->L;
	}

	vector<vector<double>> getU() {
		return this->U;
	}

	void ge_partial_pivoting(vector<vector<double>>& a, vector<double>& b) {
		this->x_GE = a;

		int n = x_GE.size();

		for (int i = 0; i < n; i++) {
			x_GE[i].push_back(b[i]); 
		}

		for (int i = 0; i < n - 1; i++) {
			int pivotRow = i;
			for (int j = i + 1; j < n; j++) {
				if (abs(x_GE[j][i]) > abs(x_GE[pivotRow][i])) {
					pivotRow = j;
				}
			}

			if (pivotRow != i) {
				swap(x_GE[i], x_GE[pivotRow]);
				swap(b[i], b[pivotRow]);
			}

			for (int j = i + 1; j < n; j++) {
				double temp = x_GE[j][i] / x_GE[i][i];
				for (int k = i; k < n + 1; k++) {
					x_GE[j][k] = x_GE[j][k] - x_GE[i][k] * temp;
				}
			}
		}

		GE.resize(n, 0.0);
		
		for (int i = n - 1; i >= 0; i--) {
			double sum = 0.0;
			for (int j = i + 1; j < n; j++) {
				sum += x_GE[i][j] * GE[j];
			}
			GE[i] = (x_GE[i][n] - sum) / x_GE[i][i];
		}
	}

	bool LU_decomp(const vector<vector<double>>& A, vector<double>& b) {
		int n = A.size();
		L.resize(n, vector<double>(n, 0.0));
		U.resize(n, vector<double>(n, 0.0));
		x_LU.resize(n, 0.0);
		y_LU.resize(n, 0.0);


		for (int i = 0; i < n; i++) {
			L[i][i] = 1.0; 
			for (int j = i; j < n; j++) {
				U[i][j] = A[i][j];
				for (int k = 0; k < i; k++) {
					U[i][j] -= L[i][k] * U[k][j];
				}
			}
			for (int j = i + 1; j < n; j++) {
				L[j][i] = A[j][i];
				for (int k = 0; k < i; k++) {
					L[j][i] -= L[j][k] * U[k][i];
				}
			
				if (U[i][i] == 0.0) {
					cout << "Chyba. Deleni nulou." << endl;
					return false;
				}
				L[j][i] /= U[i][i];
			}
		}


		for (int i = 0; i < n; i++) {
			double sum = 0.0;
			for (int j = 0; j < i; j++) {
				sum += L[i][j] * y_LU[j];
			}
			y_LU[i] = b[i] - sum;
		}


		for (int i = n - 1; i >= 0; i--) {
			double sum = 0.0;
			for (int j = i + 1; j < n; j++) {
				sum += U[i][j] * x_LU[j];
			}

			if (U[i][i] == 0.0) {
				cout << "Chyba. Deleni nulou." << endl;
				return false;
			}
			x_LU[i] = (y_LU[i] - sum) / U[i][i];
		}
		return true;
	}

	void printMatrix(vector<vector<double>> x) {
		int n = x.size();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout << x[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}

	void printVector(vector<double> x,char id) {
		int n = x.size();
		for (int i = 0; i < n; i++) {
			cout << id << "[" << i << "] = " << x[i] << endl;
		}
		cout << endl;
	}

};

vector<vector<double>> A_1 = {
{1,2,-3,4,5},
{2,1,0,3,0},
{0,2,1,2,-1},
{3,-1,0,5,2},
{2,1,2,3,-4}
};
vector<double> B_1 = { 4,9,5,3,2 };

vector<vector<double>> A_2 = {
{2,4,-2,-2},
{1,2,4,-3},
{-3,-3,8,-2},
{-1,1,6,-3}
};
vector<double> B_2 = { -4,5,7,7 };

vector<vector<double>> A_3 = {
{1,2,3,4,5},
{6,7,8,9,10},
{11,12,13,14,15},
{16,17,18,19,20},
{21,22,23,24,25}
};
vector<double> B_3 = { 1,2,3,4,5 };

vector<vector<double>> A_4 = {
{2,5,0,8},
{1,4,2,6},
{7,8,9,3},
{1,5,7,8}
};
vector<double> B_4 = { 6,5,4,3 };

int main() {
	Matrix_and_fun matrix;

	vector<vector<double>> L;
	vector<vector<double>> U;

	matrix.setA(A_1);
	matrix.setB(B_1);

	cout << "GE 1:" << endl;
	matrix.ge_partial_pivoting(A_1, B_1);
	matrix.printMatrix(matrix.getX_ge());
	matrix.printVector(matrix.getGE(), 'x');
	cout << endl;

	matrix.setA(A_2);
	matrix.setB(B_2);

	cout << "GE 2:" << endl;
	matrix.ge_partial_pivoting(A_2, B_2);
	matrix.printMatrix(matrix.getX_ge());
	matrix.printVector(matrix.getGE(), 'x');
	cout << endl;

	matrix.setA(A_3);
	matrix.setB(B_3);

	cout << "GE 3:" << endl;
	matrix.ge_partial_pivoting(A_3, B_3);
	matrix.printMatrix(matrix.getX_ge());
	matrix.printVector(matrix.getGE(), 'x');
	cout << endl;

	matrix.setA(A_4);
	matrix.setB(B_4);

	cout << "GE 4:" << endl;
	matrix.ge_partial_pivoting(A_4, B_4);
	matrix.printMatrix(matrix.getX_ge());
	matrix.printVector(matrix.getGE(), 'x');
	cout << endl;

	cout << "LU:" << endl << "------------------------------------------------" << endl;
	if (matrix.LU_decomp(A_1, B_1) == 1) {
		cout << "Matrix L 1: " << endl;
		matrix.printMatrix(matrix.getL());
		cout << "Matrix U 2: " << endl;
		matrix.printMatrix(matrix.getU());
		matrix.printVector(matrix.get_xLU(), 'x');
		matrix.printVector(matrix.get_yLU(), 'y');
	}
	cout << endl;

	cout << "LU 2: " << endl;
	if (matrix.LU_decomp(A_2, B_2) == 1){
		cout << "Matrix L2: " << endl;
		matrix.printMatrix(matrix.getL());
		cout << "Matrix U2: " << endl;
		matrix.printMatrix(matrix.getU());
		matrix.printVector(matrix.get_xLU(), 'x');
		matrix.printVector(matrix.get_yLU(), 'y');
	}
	cout << endl;
	
	cout << "LU 3: " << endl;
	if(matrix.LU_decomp(A_3, B_3) == 1){
		cout << "Matrix L3: " << endl;
		matrix.printMatrix(matrix.getL());
		cout << "Matrix U3: " << endl;
		matrix.printMatrix(matrix.getU());
		matrix.printVector(matrix.get_xLU(), 'x');
		matrix.printVector(matrix.get_yLU(), 'y');
	}
	cout << endl;

	cout << "LU 4: " << endl;
	if (matrix.LU_decomp(A_4, B_4) == 1) {
		cout << "Matrix L 4: " << endl;
		matrix.printMatrix(matrix.getL());
		cout << "Matrix U 4: " << endl;
		matrix.printMatrix(matrix.getU());
		matrix.printVector(matrix.get_xLU(), 'x');
		matrix.printVector(matrix.get_yLU(), 'y');
	}
	cout << endl;

	return 0;
}