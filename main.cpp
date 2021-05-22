/*  
 *  Least Square Approximation Chart
 *  @author  Aliaksei Korshuk
 *  @version 1.0
 *  @since   2021-03-04
 */

#include<iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <cstring>

#ifdef WIN32
    #define GNUPLOT_NAME "D:\\gnuplot\\bin\\gnuplot -persist"
#else
    #define GNUPLOT_NAME "gnuplot -persist"
#endif

using namespace std;
class SquareMatrix;
class IdentityMatrix;
class EliminationMatrix;
class PermutationMatrix;

/**
 * Matrix class realisation
 */
class Matrix {
public:
    /**
     * Default constructor
     */
    Matrix(){
        this->rows = 0;
        this->columns = 0;
        this->matrix = nullptr;
    }

    /**
     * Constructor with parameters
     *
     * @param rows size
     * @param columns size
     * @param matrix rows*columns matrix
     */
    Matrix(int rows, int columns, double** matrix){
        this->rows = rows;
        this->columns = columns;

        double** new_matrix = new double*[rows];
        for (int x = 0; x < rows; x++){
            new_matrix[x] = new double[columns];
            for (int y = 0; y < columns; y++){
                new_matrix[x][y] = matrix[x][y];
            }
        }

        this->matrix = new_matrix;
    }

    /**
     * Constructor with parameters
     *
     * @param rows size
     * @param columns size
     */
    Matrix(int rows, int columns){
        this->rows = rows;
        this->columns = columns;

        double** new_matrix = new double*[rows];
        for (int x = 0; x < rows; x++){
            new_matrix[x] = new double[columns];
            for (int y = 0; y < columns; y++){
                new_matrix[x][y] = 0;
            }
        }
        this->matrix = new_matrix;
    }

    /**
     * Overloading operator <<
     *
     * @param out Ostream
     * @param m Matrix
     * @return Ostream
     */
    friend ostream & operator << (ostream &out, const Matrix &m){
        double** matrix = m.matrix;
        for (int x = 0; x < m.rows; x++){
            out << matrix[x][0];
            for (int y = 1; y < m.columns; y++){
                out << " " << matrix[x][y];
            }
            out << endl;
        }
        return out;
    }

    /**
     * Overloading operator >>
     *
     * @param in Istream
     * @param m Matrix
     * @return Istream
     */
    friend istream & operator >> (istream &in, Matrix &m){
        in >> m.rows >> m.columns;
        double** new_matrix = new double*[m.rows];
        for (int x = 0; x < m.rows; x++){
            new_matrix[x] = new double[m.columns];
            for (int y = 0; y < m.columns; y++){
                in >> new_matrix[x][y];
            }
        }
        m.matrix = new_matrix;
        return in;
    }

    /**
     * Overloading operator +
     *
     * @param b Matrix term of addition
     * @return Resulting matrix of the operation
     */
    Matrix operator+ (const Matrix& b){

        if (rows == b.rows && columns == b.columns){
            Matrix temp(rows, columns, matrix);
            for (int x = 0; x < rows; x++){
                for (int y = 0; y < columns; y++){
                    temp.matrix[x][y]+=b.matrix[x][y];
                }
            }
            return temp;
        }
        Matrix temp;
        cout << "Error: the dimensional problem occurred" << endl;
        return temp;
    }

    /**
     * Overloading operator -
     *
     * @param b Subtrahend Matrix
     * @return Resulting matrix of the operation
     */
    Matrix operator- (const Matrix& b){

        if (rows == b.rows && columns == b.columns){
            Matrix temp(rows, columns, matrix);
            for (int x = 0; x < rows; x++){
                for (int y = 0; y < columns; y++){
                    temp.matrix[x][y]-=b.matrix[x][y];
                }
            }
            return temp;
        }
        Matrix temp;
        cout << "Error: the dimensional problem occurred" << endl;
        return temp;
    }
    /**
     * Overloading operator *
     *
     * @param b Factor Matrix
     * @return Resulting matrix of the operation
     */
    Matrix operator* (const Matrix& b){

        if (columns == b.rows){
            Matrix temp(rows, b.columns);
            for (int x = 0; x < rows; x++){
                for (int y = 0; y < b.columns; y++){
                    for(int z=0; z<columns; z++) {
                        temp.matrix[x][y]+=matrix[x][z]*b.matrix[z][y];
                    }
                }
            }
            return temp;
        }
        Matrix temp;
        cout << "Error: the dimensional problem occurred" << endl;
        return temp;
    }

    /**
     * Calculating transpose matrix
     *
     * @return Resulting matrix of the operation
     */
    Matrix transpose()
    {
        Matrix temp(columns, rows);
        for (int i = 0; i < columns; i++)
            for (int j = 0; j < rows; j++)
                temp.matrix[i][j] = matrix[j][i];
        return temp;
    }

public:
    int rows;       // Amount of rows
    int columns;    // Amount of columns
    double** matrix;   // Matrix of size rows*columns
};

/**
 * SquareMatrix class realisation
 */
class SquareMatrix:public Matrix {
public:
    /**
     * Default constructor
     */
    SquareMatrix() : Matrix(){
    }

    /**
     * Constructor with parameters
     *
     * @param rows size
     * @param matrix matrix size*size
     */
    SquareMatrix(int size, double** matrix) : Matrix(size,size, matrix){
    }

    /**
     * Constructor with parameters
     *
     * @param size size
     */
    SquareMatrix(int size): Matrix(){
        this->rows = size;
        this->columns = size;

        double** new_matrix = new double*[size];
        for (int x = 0; x < size; x++){
            new_matrix[x] = new double[size];
            for (int y = 0; y < size; y++){
                new_matrix[x][y] = 0;
            }
        }
        this->matrix = new_matrix;
    }

    /**
     * Overloading operator <<
     *
     * @param out Ostream
     * @param m Matrix
     * @return Ostream
     */
    friend ostream & operator << (ostream &out, const SquareMatrix &m){
        double** matrix = m.matrix;
        for (int x = 0; x < m.rows; x++){
            out << matrix[x][0];
            for (int y = 1; y < m.rows; y++){
                out << " " << matrix[x][y];
            }
            out << endl;
        }
        return out;
    }

    /**
     * Overloading operator >>
     *
     * @param in Istream
     * @param m Matrix
     * @return Istream
     */
    friend istream & operator >> (istream &in, SquareMatrix &m){
        in >> m.rows;
        m.columns = m.rows;
        double** new_matrix = new double*[m.rows];
        for (int x = 0; x < m.rows; x++){
            new_matrix[x] = new double[m.rows];
            for (int y = 0; y < m.rows; y++){
                in >> new_matrix[x][y];
            }
        }
        m.matrix = new_matrix;
        return in;
    }

    /**
     * Overloading operator +
     *
     * @param b Matrix term of addition
     * @return Resulting matrix of the operation
     */

    SquareMatrix operator+ (const SquareMatrix& b){
        Matrix newMatrix = (Matrix)*this + (Matrix)b;
        return *(SquareMatrix*)(&newMatrix);
    }

    /**
     * Overloading operator -
     *
     * @param b Subtrahend Matrix
     * @return Resulting matrix of the operation
     */
    SquareMatrix operator- (const SquareMatrix& b){
        Matrix newMatrix = (Matrix)*this - (Matrix)b;
        return *(SquareMatrix*)(&newMatrix);
    }
    /**
     * Overloading operator *
     *
     * @param b Factor Matrix
     * @return Resulting matrix of the operation
     */
    SquareMatrix operator* (const SquareMatrix& b){
        Matrix newMatrix = (Matrix)*this * (Matrix)b;
        return *(SquareMatrix*)(&newMatrix);
    }

    /**
     * Calculating transpose matrix
     *
     * @return Resulting matrix of the operation
     */
    SquareMatrix transpose()
    {
        Matrix newMatrix = (Matrix)*this;
        newMatrix = newMatrix.transpose();
        return *(SquareMatrix*)(&newMatrix);
    }

};


/**
 * IdentityMatrix class realisation
 */
class IdentityMatrix : public SquareMatrix{
public:

    /**
     * Constructor with parameters
     *
     * @param size Size of matrix
     */
    IdentityMatrix(int size): SquareMatrix(size){
        for (int x = 0; x < size; x++){
            this->matrix[x][x] = 1;
        }
    }

};

/**
 * EliminationMatrix class realisation
 */
class EliminationMatrix  : public IdentityMatrix{
public:

    /**
     * Constructor with parameters
     *
     * @param x Position
     * @param y Position
     * @param matrix Input matrix
     */
    EliminationMatrix(int x, int y, SquareMatrix matrix): IdentityMatrix(matrix.rows){
        double value;
        double a = matrix.matrix[x - 1][y - 1];
        double k = matrix.matrix[y - 1][y - 1];

        value = -1 * a / k;

        this->matrix[x - 1][y - 1] = value;
    }
};

/**
 * PermutationMatrix class realisation
 */
class PermutationMatrix  : public IdentityMatrix{
public:

    /**
     * Constructor with parameters
     *
     * @param x Position
     * @param y Position
     * @param matrix Input matrix
     */
    PermutationMatrix(int x, int y, SquareMatrix matrix): IdentityMatrix(matrix.rows){
        double* temp;

        temp = this->matrix[x-1];
        this->matrix[x-1] = this->matrix[y - 1];
        this->matrix[y - 1] = temp;
    }
};



/**
 * Method that swaps rows by index
 * @param x Index
 * @param y Index
 * @param A SquareMatrix
 */
void swap( int x, int y, SquareMatrix A)
{
    double* temp;
    temp = A.matrix[x];
    A.matrix[x] = A.matrix[y];
    A.matrix[y] = temp;
}

/**
 * Calculates the determinant after Gaussian elimination
 *
 * @param A SquareMatrix
 * @return determinant
 */
double determinant(SquareMatrix A){
    double answer = 1;
    for (int i = 0; i < A.rows; i++) answer*=A.matrix[i][i];
    return answer;
}

/**
 * Method that calculates Inverse Matrix using Gaussian elimination
 * @param A SquareMatrix
 */
SquareMatrix inverseMatrix(SquareMatrix A){
    SquareMatrix I = IdentityMatrix(A.rows);

    int count = 1;
    for (int k=0; k<A.rows; k++)
    {
        int i_max = k;
        double v_max = abs(A.matrix[i_max][k]);

        for (int i = A.rows-1; i > k; i--)
            if (abs(A.matrix[i][k]) > v_max)
                v_max = A.matrix[i][k], i_max = i;


        if (i_max != k){
            swap(k, i_max, A);
            swap(k, i_max, I);

            count++;
        }


        for (int i=k+1; i<A.rows; i++)
        {

            if (A.matrix[i][k] != 0){
                EliminationMatrix E(i + 1, k + 1, A);
                A = E*A;
                I = E*I;

                count++;
            }

        }
    }

    for (int k=A.rows-1; k>0; k--)
    {
        for (int i=A.rows-2 - (A.rows - 1 - k); i>=0; i--)
        {

            if (A.matrix[i][k] != 0){
                EliminationMatrix E(i + 1, k + 1, A);
                A = E*A;
                I = E*I;
                count++;
            }

        }
    }

    for (int i = 0; i < A.rows; i++){
        double temp = 1/A.matrix[i][i];
        for (int k = 0; k < I.rows; k++){
            I.matrix[i][k]*=temp;
        }
        A.matrix[i][i] = 1;
    }

    return I;

}


/**
 * Method that solves Least Square Approximation
 *
 * @param A Matrix
 * @param b Matrix
 */
Matrix leastSquareApproximation(Matrix A, Matrix b){
    cout << "A:\n" << A;
    cout << "A_T*A:\n" << A.transpose()*A;
    Matrix temp = A.transpose()*A;
    Matrix a_new = inverseMatrix(*(SquareMatrix*)(&temp));
    cout << "(A_T*A)^-1:\n" << a_new;
    Matrix b_new = A.transpose()*b;
    cout << "A_T*b:\n" << b_new;
    cout << "x~:\n";

    return a_new*b_new;

}

// Driver Code
int main() {
    cout << fixed << setprecision(4);
    srand(time(0));
    int size;
    size = 20;
    double* t = new double[size];
    Matrix b(size, 1);

    for (int i = 0; i < size; i++){
        t[i] = i+1;
        b.matrix[i][0] = rand() % 4 + i;
    }

    int degree;
    degree = 10;
    Matrix A(size, degree + 1);
    for (int i = 0; i < size; i++){
        for (int k = 0; k <= degree; k++){
            A.matrix[i][k] = pow(t[i], k);
        }
    }

    Matrix X = leastSquareApproximation(A, b);
    cout << X;

    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    string func;
    for (int i = 0; i <= degree; i++){
        stringstream ss;
        ss << X.matrix[i][0] ;
        string s;
        ss >> s;
        func+= s;
        func+= "*x**";
        stringstream ss2;
        ss2 << i ;
        string s2;
        ss2 >> s2;
        func+= s2;
        func+=" + ";
    }
    func+="0";
    cout << func << endl;

    int n = func.length();

    char char_array[n + 1];

    strcpy(char_array, func.c_str());

    fprintf(pipe, "set yrange [0:%d]\n", size + 6);
    fprintf(pipe, "f(x) = %s\n", char_array);

    fprintf(pipe, "plot [0:%d] f(x) title 'appr', '-' using 1:2 title 'exp' with points pointtype 5 pointsize 1\n", size);
    for (int i = 0; i < size; i++){
        fprintf(pipe, "%f %f\n", t[i], b.matrix[i][0]);
    }

    fprintf(pipe, "%s\n", "e");

    _pclose(pipe);

    return 0;
}
