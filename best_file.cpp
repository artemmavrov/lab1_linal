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
	void write_matrix()
	{
		cin >> this->n >> this->m;
		rowpointer.resize(n + 1, 0);
		for (int i = 0; i < this->n; i++)
		{
			int cnt_of_not_zero_elements = 0;
			for (int j = 0; j < this->m; j++)
			{
				long double value;
				cin >> value;
				if (value != 0)
				{
					this->values.push_back(value);
					this->col.push_back(j);
					cnt_of_not_zero_elements++;
				}
			}
			rowpointer[i + 1] = rowpointer[i] + cnt_of_not_zero_elements;
		}
	}

public:
	int n, m;
	long double scalar = 1.0;
	vector<long double> values;
	vector<int> col;
	vector<int> rowpointer;
	string path;

	Matrix() = default;
	Matrix(int n = 0, int m = 0, vector<long double> values = {}, vector<int> col = {}, vector<int> rowpointer = {}, string path = "") {
        this->n = n;
        this->m = m;
        this->path = path;

        if (!values.empty()) {
            this->values = values;
            this->col = col;
            this->rowpointer = rowpointer;
        } else if (!path.empty()) {
			
            ifstream file(path);
			
            file >> this->n >> this->m;
            long double value;
			this->rowpointer.resize(this->n + 1, 0);
            for (int i = 0; i < this->n; ++i) {
				int cnt_of_not_zero_elements = 0;
                for (int j = 0; j < this->m; ++j) {
                    file >> value;
                    if (value != 0) {
                        this->values.push_back(value);
                        this->col.push_back(j);
						cnt_of_not_zero_elements++;
                    }
                }
				this->rowpointer[i + 1] = this->rowpointer[i] + cnt_of_not_zero_elements;
            }
            file.close();
			
        } else {
            write_matrix();
        }
    }

	long double get_trace() {
		long double trace = 0;
		for (int i = 1; i < rowpointer.size(); i++)
		{
			for(int j = rowpointer[i - 1]; j < rowpointer[i]; j++)
			{
				// rowpointer[i] - right border, rowpointer[i - 1] - left border
				// if the diagonal element by j coord not in [rowpointer[i - 1], rowpointer[i]) we can skip action
				if(rowpointer[i] <= i || rowpointer[i - 1] > i)
				{
					continue;
				}
				if(col[j] == i - 1)
				{
					trace += values[j];
				}
			}
		}

		return trace * scalar;
	}

	long double get_element(int i, int j)
	{
		j--; // indexation, don't write i--;!!!
		// if element not in our borders, it means that element value = 0
		if(col[rowpointer[i - 1]] > j || col[rowpointer[i] - 1] < j)
		{
			return 0;
		}
		for(int value_index = rowpointer[i - 1]; value_index < rowpointer[i]; value_index++)
		{
			if(col[value_index] == j)
			{
				return values[value_index];
			}
		}
		return 0;
    }

	//function that transforms a sparse matrix into a dense one
	vector<vector<long double>> to_vector() const {
        vector<vector<long double>> dense_matrix(n, vector<long double>(m, 0));

		for(int i = 1; i < rowpointer.size(); i++)
		{
			for(int j = rowpointer[i - 1]; j < rowpointer[i]; j++)
			{
				dense_matrix[i - 1][col[j]] = values[j] * scalar;
			}
		}

        return dense_matrix;
    }
};

// TODO: rename this function
void get_matrix_column_sum(map<int, long double> &column_values, Matrix &matrix, int &col_index)
{
	for(int row_index = 0; row_index < matrix.n; row_index++)
	{
		long double element = matrix.get_element(row_index + 1, col_index + 1);
		if(element != 0)
		{
			column_values[row_index] = element;
		}
	}
}

// TODO: rename this function
void get_matrix_string_sum(map<int, long double> &string_values, Matrix &matrix, int &rowpointer_index)
{
		for(int values_index = matrix.rowpointer[rowpointer_index - 1]; values_index < matrix.rowpointer[rowpointer_index]; values_index++)
		{
			string_values[matrix.col[values_index]] += matrix.values[values_index];
		}
}

void merge_string_to_new_vectors(vector<long double> &new_values, vector<int> &new_col, map<int, long double> &string_values)
{
	for(pair<int, long double> string_elements : string_values)
	{
		int col = string_elements.first;
		long double value = string_elements.second;
		new_values.push_back(value);
		new_col.push_back(col);
	}
}

Matrix sum_of_matrix(Matrix& A, Matrix& B) {
	if (A.n != B.n || A.m != B.m) {
		cout << "Mistake";
		return A;
	}
	int n = A.n;
	int m = A.m;
	vector<long double> new_values;
	vector<int> new_col;
	vector<int> new_rowpointer(n + 1, 0);
	for(int rowpointer_index = 1; rowpointer_index < n + 1; rowpointer_index++)
	{
		map<int, long double> string_values;
		get_matrix_string_sum(string_values, A, rowpointer_index);
		get_matrix_string_sum(string_values, B, rowpointer_index);
		merge_string_to_new_vectors(new_values, new_col, string_values);
		new_rowpointer[rowpointer_index] = new_rowpointer[rowpointer_index - 1] + string_values.size();
	}
	Matrix C = Matrix(n, m, new_values, new_col, new_rowpointer);
	return C;
}

void multiply_scalar(int scalar, Matrix& A) {
	A.scalar *= scalar;
}

Matrix multiply_matrix(Matrix& A, Matrix& B) {
	if (A.m != B.n) {
		cout << "Mistake \n";
		return A;
	}
	int new_n = A.n;
	int new_m = B.m;
	vector<long double> new_values;
	vector<int> new_col;
	vector<int> new_rowpointer(new_n + 1, 0);
	// column values of matrix B
	map<int, map<int, long double>> column_values;
	// string values of matrix A
	map<int, map<int, long double>> string_values;
	for(int i = 1; i < A.rowpointer.size(); i++)
	{
		// TODO: check this code later
		get_matrix_string_sum(string_values[i - 1], A, i);
	}
	for(int j = 0; j < B.m; j++)
	{
		get_matrix_column_sum(column_values[j], B, j);
	}
	for(int i = 0; i < new_n; i++)
	{
		int cnt_of_not_zero_elements = 0;
		for(int j = 0; j < new_m; j++)
		{
			long double new_value = 0;
			for(pair<int, long double> string_value : string_values[i])
			{
				long double value = string_value.second;
				long double index = string_value.first;
				if(column_values[j].count(index) == 0)
				{
					continue;
				}
				value *= column_values[j][index];
				new_value += value;
			}
			new_value = new_value * A.scalar * B.scalar;
			
			if(new_value != 0)
			{
				cnt_of_not_zero_elements++;
				new_values.push_back(new_value);
				new_col.push_back(j);
			}
		}
		new_rowpointer[i + 1] = new_rowpointer[i] + cnt_of_not_zero_elements;
	}
	Matrix C = Matrix(new_n, new_m, new_values, new_col, new_rowpointer);
	return C;
}

pair<long double, string> get_determinant(Matrix matrix) {
	int n = matrix.n;
	vector<long double> values = matrix.values;
	vector<int> col = matrix.col;
	vector<int> rowpointer = matrix.rowpointer;

	vector<int> order(n + 1);
	vector<long double> nw_values;
	vector<int> nw_col;
	vector<int> nw_rowpointer;
	long double det = 1.0;
	long double factor;
	int pivot_row;
	long double max_value;

	for (int i = 0; i < n + 1; i++) {
		order[i] = i;
	}

	for (int i = 0; i < n; i++) {
		max_value = 0.0;

		for (int j = i; j < n; j++) {
			for (int k = rowpointer[order[j]]; k < rowpointer[order[j] + 1]; k++) {
				if (col[k] == i) {
					if (abs(values[k]) > abs(max_value)) {
						max_value = values[k];
						pivot_row = j;
					}

					break;
				}
			}
		}

		if (max_value == 0.0) {
			return make_pair(0.0, "No");
		}

		det *= max_value;

		if (pivot_row != i) {
			det *= -1;
			swap(order[pivot_row], order[i]);						
		}

		for (int j = i + 1; j < n; j++) {
			if (col[rowpointer[order[j]]] == i) {
				factor = values[rowpointer[order[j]]] / max_value;

				values[rowpointer[order[j]]] = 0.0;
				for (int k = rowpointer[order[j]] + 1, l = rowpointer[order[i]] + 1; k < rowpointer[order[j] + 1]; k++) {
					if (col[k] == col[l]) {
						values[k] -= values[l] * factor;
						l++;
					}
				}
			}
		}

		nw_rowpointer = rowpointer;
		for (int j = 0; j < n; j++) {
			for (int k = rowpointer[j]; k < rowpointer[j + 1]; k++) {
				if (abs(values[k]) > 1e-5) {
					nw_values.push_back(values[k]);
					nw_col.push_back(col[k]);
				}else {
					for (int p = j + 1; p < n + 1; p++) {
						nw_rowpointer[p] -= 1;
					}
				}
			}
		}
		values = move(nw_values);
		col = move(nw_col);
		rowpointer = move(nw_rowpointer);
	}

	return make_pair(det, "Yes");
}

