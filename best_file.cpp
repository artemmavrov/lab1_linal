#include <iostream>
#include <map>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

class Matrix
{
private:
	// ручной ввод матрицы
	void write_matrix()
	{
		cin >> this->n >> this->m;
		for (int i = 0; i < this->n; i++)
		{
			for (int j = 0; j < this->m; j++)
			{
				long double element;
				cin >> element;
				if (element != 0)
				{
					this->elements.push_back(element);
					this->indexI.push_back(i);
					this->indexJ.push_back(j);
				}
			}
		}
	}

public:
	int n, m;
	long double scalar = 1.0;
	vector<long double> elements;
	vector<int> indexI;
	vector<int> indexJ;
	string path;

	Matrix() = default;
	Matrix(int n = 0, int m = 0, vector<long double> elements = {}, vector<int> indexI = {}, vector<int> indexJ = {}, string path = "") {
        this->n = n;
        this->m = m;
        this->path = path;

        //input either via file or via console (if both options are submitted, file is a priority)
        if (!elements.empty()) {
            this->elements = elements;
            this->indexI = indexI;
            this->indexJ = indexJ;
        } else if (!path.empty()) {
            ifstream file(path);

            file >> this->n >> this->m;
            long double element;

            for (int i = 0; i < this->n; ++i) {
                for (int j = 0; j < this->m; ++j) {
                    file >> element;
                    if (element != 0) {
                        this->elements.push_back(element);
                        this->indexI.push_back(i);
                        this->indexJ.push_back(j);
                    }
                }
            }

            file.close();
        } else {
            write_matrix();
        }
    }

	long double get_trace() {
		long double trace = 0;
		for (int i = 0; i < elements.size(); i++) {
			if (indexI[i] == indexJ[i]) {
				trace += elements[i];
			}
		}

		return trace * scalar;
	}

	long double get_element(int i, int j) {
      i--; j--; // indexation
      for (int idx = 0; idx < elements.size(); idx++)
      {
        if (indexI[idx] == i && indexJ[idx] == j)
        {
          return elements[idx] * scalar;
        }
      }
      return 0;
    }

	//function that transforms a sparse matrix into a dense one
	vector<vector<long double>> to_vector() const {
        vector<vector<long double>> dense_matrix(n, vector<long double>(m, 0));

        for (size_t i = 0; i < elements.size(); ++i) {
            dense_matrix[indexI[i]][indexJ[i]] = elements[i] * scalar;
        }

        return dense_matrix;
    }
};

// вспомогательная функция для sum_of_matrix
void write_map(Matrix& A, map<pair<int, int>, long double>& map) {
	for (int i = 0; i < A.elements.size(); i++)
	{
		pair<int, int> coord = { A.indexI[i], A.indexJ[i] };
		map[coord] += A.elements[i] * A.scalar;
	}
}

Matrix sum_of_matrix(Matrix& A, Matrix& B) {
	if (A.n != B.n || A.m != B.m) {
		// Обработать описание ошибки
		cout << "Mistake";
		return A;
	}
	int n = A.n;
	int m = A.m;
	vector<long double> elements;
	vector<int> indexI;
	vector<int> indexJ;
	int indA = 0; // indicator A - указатель
	int indB = 0; // indicator B - указатель
	int end = max(A.elements.size(), B.elements.size());

	// С помощью unordered_map
	map<pair<int, int>, long double> map;
	write_map(A, map);
	write_map(B, map);

	for (auto map_element : map)
	{
		pair<int, int> coord = map_element.first;
		int matrix_element = map_element.second;
		elements.push_back(matrix_element);
		indexI.push_back(coord.first);
		indexJ.push_back(coord.second);
	}

	Matrix C = Matrix(n, m, elements, indexI, indexJ);
	return C;
}

void multiply_scalar(int scalar, Matrix& A) {
	A.scalar *= scalar;
}

vector<long double> get_vertical(Matrix& A, int j) {
	vector<long double> vertical(A.n, 0);

	for (int idx = 0; idx < A.elements.size(); idx++)
	{
		// Если нашли элемент по горизонтали, то в соответствующую вертикальную ячейку записываем сам элемент
		if (A.indexJ[idx] == j)
		{
			// т.к. вектор вертикальный, то и обращаться нужно к A.indexI, т.к. этот вектор отвечает за вертикаль
			vertical[A.indexI[idx]] = A.elements[idx] * A.scalar;
		}
	}

	return vertical;
}

vector<long double> get_horizontal(Matrix& A, int i) {
	vector<long double> horizontal(A.m, 0);

	for (int idx = 0; idx < A.elements.size(); idx++) {
		// Если нашли элемент по вертикали, то в соответствующую горизонтальную ячейку записываем сам элемент
		if (A.indexI[idx] == i) {
			// т.к. вектор горизонтальный, то и обращаться нужно к A.indexJ, т.к. этот вектор отвечает за горизонталь
			horizontal[A.indexJ[idx]] = A.elements[idx] * A.scalar;
		}
	}

	return horizontal;
}

long double multiply_vectors(vector<long double>& a, vector<long double>& b) {
	long double res = 0;

	for (long double i = 0; i < a.size(); i++) {
		res += a[i] * b[i];
	}

	return res;
}

Matrix multiply_matrix(Matrix& A, Matrix& B) {
	if (A.m != B.n) {
		// Обработать описание ошибки
		cout << "Mistake \n";
		return A;
	}

	int new_n = A.n;
	int new_m = B.m;
	vector<long double> elements;
	vector<int> indexI;
	vector<int> indexJ;

	for (int i = 0; i < new_n; i++)
	{
		// что бы не пересчитывать каждый раз, вынес horizontal
		vector<long double> horizontal = get_horizontal(A, i);
		for (int j = 0; j < new_m; j++)
		{
			vector<long double> vertical = get_vertical(B, j);
			long double composition = multiply_vectors(horizontal, vertical);
			if (composition != 0)
			{
				elements.push_back(composition);
				indexI.push_back(i);
				indexJ.push_back(j);
			}
		}
	}

	Matrix C = Matrix(new_n, new_m, elements, indexI, indexJ);
	return C;
}

pair<long double, string> get_determinant(Matrix matrix) {
    vector<vector<long double>> temp_matrix = matrix.to_vector();
    int rows = matrix.n;

    long double det = 1.0;

    for (int i = 0; i < rows; ++i) {
        int pivot_row = i;

        for (int j = i + 1; j < rows; ++j) {
            if (abs(temp_matrix[j][i]) > abs(temp_matrix[pivot_row][i])) {
                pivot_row = j;
            }
        }

        if (abs(temp_matrix[pivot_row][i]) < 1e-9) {
            return make_pair(0.0, "No");
        }

        if (pivot_row != i) {
            swap(temp_matrix[i], temp_matrix[pivot_row]);
            det = -det; // Меняем знак детерминанта при перестановке строк
        }

        det *= temp_matrix[i][i];

        for (int j = i + 1; j < rows; ++j) {
            long double factor = temp_matrix[j][i] / temp_matrix[i][i];
            for (int k = i; k < rows; ++k) {
                temp_matrix[j][k] -= factor * temp_matrix[i][k];
            }
        }
    }

    if (det != 0.0) {
        return make_pair(det, "Yes");
    } else {
        return make_pair(det, "No");
    }
}
