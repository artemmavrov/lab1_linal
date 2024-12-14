#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <windows.h>
#include <cmath>
#include "best_file.cpp"

using namespace std;

class dense_matrix {
private:
    int rows, cols;
    vector<vector<long double>> matrix;

public:
    dense_matrix(const string& filename) {
        ifstream file(filename);

        file >> rows >> cols;
        matrix.resize(rows, vector<long double>(cols));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                file >> matrix[i][j];
            }
        }

        file.close();
    }

    long double get_tace_dense() const {
        if (rows != cols) {
            throw invalid_argument("Trace is defined only for square matrices");
        }

        int trace = 0;
        for (int i = 0; i < rows; ++i) {
            trace += matrix[i][i];
        }

        return trace;
    }

    long double get_element_dense(int i, int j) const {
        return matrix[i - 1][j - 1];
    }

    dense_matrix addition_matrix_dense(const dense_matrix& other) const {
        dense_matrix result(rows, cols);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
            }
        }

        return result;
    }

    dense_matrix multiply_scalar_dense(long double scalar) const {
        dense_matrix result(rows, cols);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result.matrix[i][j] = matrix[i][j] * scalar;
            }
        }

        return result;
    }

    dense_matrix multiply_matrix_dense(const dense_matrix& other) const {
        dense_matrix result(rows, other.cols);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                result.matrix[i][j] = 0;
                for (int k = 0; k < cols; ++k) {
                    result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
                }
            }
        }

        return result;
    }

    dense_matrix(int r, int c) : rows(r), cols(c) {
        matrix.resize(rows, vector<long double>(cols, 0));
    }

     vector<vector<long double>> to_vector_dense() const {
        return matrix;
    }

    pair<long double, string> get_determinant_dense() const {
        vector<vector<long double>> temp_matrix = matrix;
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
                det = -det;
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
};

string buffer;

void write_matrix(vector<vector<long double>>& matrix, int& n, int& m) {
	for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
	}
}

void generate_matrix(int& n, int& m, int& perecent_null, int number_matrix) {
    auto seed = chrono::system_clock::now().time_since_epoch().count() + number_matrix;
    mt19937 gen(seed);
    uniform_real_distribution<long double> rand_element(-1e3, 1e3);
    uniform_int_distribution<int> rand_perecent(0, 100);

    ofstream test;
    string file_way = buffer + "/tests/matrix_" + to_string(number_matrix) + ".txt";
    test.open(file_way);

    test << n << " " <<  m << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (rand_perecent(gen) < perecent_null) {
                test << 0 << " ";
            }else {
                test << rand_element(gen) << " ";
            }
        }
        test << endl;
    }

    test.close();
}

bool check_answer_double(long double& answer, long double& correct_answer, long double& error) {
    return fabs(abs(answer - correct_answer)) <= error;
}

bool check_answer_matrix(vector<vector<long double>>& answer, vector<vector<long double>>& correct_answer, int& n, int& m, long double& error) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (!check_answer_double(answer[i][j], correct_answer[i][j], error)) {
                return false;
            }
        }
    }

    return true;
}

void run_test(int& n, int& m, int& number_martrix, int& testing_func, long double& error) {
    string file_way_matrix1 = buffer + "/tests/matrix_" + to_string(number_martrix) + ".txt";
    Matrix matrix1 = Matrix(0, 0, {}, {}, {}, file_way_matrix1);
    dense_matrix matrix1_dense = dense_matrix(file_way_matrix1);

    if (testing_func != 3 && testing_func != 5) {
      cout << "Test " + to_string(number_martrix + 1) + ": ";
    }else {
        cout << "Test " + to_string(number_martrix / 2 + 1) + ": ";
    }

    if (testing_func == 1) {
        long double answer = matrix1.get_trace();
        long double correct_answer = matrix1_dense.get_tace_dense();

         if (check_answer_double(answer, correct_answer, error)) {
            cout << "Accept" << endl;
         }else {
            cout << "Wrong answer"
            << endl << " Your code output: " << answer
            << endl << " Correct answer: " << correct_answer << endl;
         }
    }else if (testing_func == 2) {
        bool check_flag = true;
        long double answer, correct_answer;

        for (int i = 1 ; i <= n && check_flag; i++) {
            for (int j = 1; j <= m; j++) {
                answer = matrix1.get_element(i, j);
                correct_answer = matrix1_dense.get_element_dense(i, j);
                if (!check_answer_double(answer, correct_answer, error)) {
                    cout << "Wrong answer on index " << i << " " << j
                    << endl << " Your code output: " << matrix1.get_element(i, j)
                    << endl << " Correct answer: " << matrix1_dense.get_element_dense(i, j) << endl;
                    check_flag = false;
                    break;
                }
            }
        }

        if (check_flag) {
            cout << "Accept" << endl;
         }
    }else if (testing_func == 3) {
        string file_way_matrix2 = buffer + "/tests/matrix_" + to_string(number_martrix + 1) + ".txt";
        Matrix matrix2 = Matrix(0, 0, {}, {}, {}, file_way_matrix2);
        dense_matrix matrix2_dense = dense_matrix(file_way_matrix2);

        vector<vector<long double>> answer = sum_of_matrix(matrix1, matrix2).to_vector();
        vector<vector<long double>> correct_answer = (matrix1_dense.addition_matrix_dense(matrix2_dense)).to_vector_dense();

        if (check_answer_matrix(answer, correct_answer, n, m, error)) {
            cout << "Accept" << endl;
         }else {
            cout << "Wrong answer" << endl << " Your code output: ";
            write_matrix(answer, n, m);
            cout <<  endl << " Correct answer: ";
            write_matrix(correct_answer, n, m);
         }
    }else if (testing_func == 4) {
        auto seed = chrono::system_clock::now().time_since_epoch().count();
        mt19937 gen(seed);
        uniform_real_distribution<long double> rand_element(-1e6, 1e6);
        long double rd_scalar = rand_element(gen);

        multiply_scalar(rd_scalar, matrix1);
        vector<vector<long double>> answer = matrix1.to_vector();
        vector<vector<long double>> correct_answer = (matrix1_dense.multiply_scalar_dense(rd_scalar)).to_vector_dense();

        if (check_answer_matrix(answer, correct_answer, n, m, error)) {
            cout << "Accept" << endl;
         }else {
            cout << "Wrong answer" << endl << " Your code output: ";
            write_matrix(answer, n, m);
            cout <<  endl << " Correct answer: ";
            write_matrix(correct_answer, n, m);
         }
    }else if (testing_func == 5) {
        string file_way_matrix2 = buffer + "/tests/matrix_" + to_string(number_martrix + 1) + ".txt";
        Matrix matrix2 = Matrix(0, 0, {}, {}, {}, file_way_matrix2);
        dense_matrix matrix2_dense = dense_matrix(file_way_matrix2);

        vector<vector<long double>> answer = multiply_matrix(matrix1, matrix2).to_vector();
        vector<vector<long double>> correct_answer = (matrix1_dense.multiply_matrix_dense(matrix2_dense)).to_vector_dense();

        if (check_answer_matrix(answer, correct_answer, n, m, error)) {
            cout << "Accept" << endl;
         }else {
            cout << "Wrong answer" << endl << " Your code output: ";
            write_matrix(answer, n, m);
            cout <<  endl << " Correct answer: ";
            write_matrix(correct_answer, n, m);
         }
    }else {
        pair<long double, string> answer = get_determinant(matrix1);
        pair<long double, string> correct_answer = matrix1_dense.get_determinant_dense();

         if (check_answer_double(answer.first, correct_answer.first, error)) {
            cout << "Accept" << endl;
         }else {
            cout << "Wrong answer"
            << endl << " Your code output: " << answer.first << " " << answer.second
            << endl << " Correct answer: " << correct_answer.first << " " << correct_answer.second << endl;
         }
    }
}

void testing() {
    int testing_function;
    cout << "Testing function(cin number):"
    << endl << "1) Matrix trace"
    << endl << "2) Element of matrix"
    << endl << "3) Matrix addition"
    << endl << "4) Matrix multiplication by scalar"
    << endl << "5) Matrix multiplication"
    << endl << "6) Matrix determinant" << endl;
    cin >> testing_function;

    int quantity;
    cout << "Quantity of tests: " << endl;
    cin >> quantity;

    int n, m;
    cout << "Height, weight: " << endl;
    cin >> n >> m;

    int perecent_null;
    cout << "Percent of null: ";
    cin >> perecent_null;
    cout << endl;

    long double error;
    cout << "Permissible error: ";
    cin >> error;
    cout << endl;

    if (testing_function == 3 || testing_function == 5) {
        for (int number_matrix = 0; number_matrix < quantity * 2; number_matrix += 2) {
            generate_matrix(n, m, perecent_null, number_matrix);
            generate_matrix(n, m, perecent_null, number_matrix + 1);
            run_test(n, m, number_matrix, testing_function, error);
        }
    }else {
        for (int number_matrix = 0; number_matrix < quantity; number_matrix++) {
            generate_matrix(n, m, perecent_null, number_matrix);
            run_test(n, m, number_matrix, testing_function, error);
        }
    }
}

int main() {
    char temp_buffer[MAX_PATH];
    GetCurrentDirectory(MAX_PATH, temp_buffer);
    buffer = string(temp_buffer);

    testing();
}

