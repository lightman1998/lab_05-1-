#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sstream>

using namespace std;

class Strategy {

    vector<vector<double>> Matrix;
    vector<vector<double >> MatrixB;
public:
    const vector<vector<double>> &getMatrixForB() const;

private:
    int size_x,size_y;
    vector<double> vector_x, vector_y;
    double min_x, max_y;
    bool saddle_point = false;


public:
    Strategy (vector<vector<double>> M){
        Matrix = M;
        size_y = M.size();
        size_x = M.at(0).size();
        minimax_x();
        maximin_y();
        isSaddlePoint();
        createMatrixForB();
    }

    double getMinX() const;

    double getMaxY() const;


    bool isSaddlePoint() {
        if(max_y == min_x){
            saddle_point = true;
            return saddle_point;
        }
        else {
            saddle_point = false;
        }
        return saddle_point;
    }

    // СТРАТЕГИЯ ИГРОКА А
    void maximin_y(){

        for (int i = 0; i < size_y; ++i) {
            vector<double> v = Matrix.at(i);
            auto result = min_element(v.begin(), v.end());
            double a = v.at(std::distance(v.begin(), result));
            vector_y.emplace_back(a);
        }

        auto max_result = max_element(vector_y.begin(), vector_y.end());
        double a = vector_y.at(std::distance(vector_y.begin(), max_result));
        max_y = a;
    }


    //СТРАТЕГИЯ ИГРОКА B
    void minimax_x(){

        for (int j = 0; j < size_x; ++j) {
            vector<double> v;
            for (int i = 0; i < size_y; ++i) {
                v.emplace_back(Matrix.at(i).at(j));
            }
            auto result = max_element(v.begin(), v.end());
            double a = v.at(std::distance(v.begin(), result));
            vector_x.emplace_back(a);
        }

        vector<double> v = vector_x;

        auto min_result = std::min_element(v.begin(), v.end());
        double a = v.at(std::distance(v.begin(), min_result));
        min_x = a;
    }


    void print();

    const vector<vector<double>> &getMatrixForA() const;
    vector<vector<double>> &createMatrixForB();


};

void Strategy::print() {

    for (int i = 0; i < size_y; ++i) {
        cout << "+";
        for (int j = 0; j < size_x + 1; ++j) {
            cout << "-------+";
        }
        cout << endl;
        for (int j = 0; j < size_x; ++j) {
            std::cout << "|  " << Matrix.at(i).at(j) << "\t";
        }
        cout << "|  " << vector_y.at(i) << "\t|\n";
    }
    cout << "+";
    for (int j = 0; j < size_x + 1; ++j) {

        cout << "-------+";
    }
    cout << endl;
    for (int k = 0; k < size_x; ++k) {
        std::cout << "|  " << vector_x.at(k) << "\t";
    }
    cout << "|\n";

    cout << "+";
    for (int j = 0; j < size_x ; ++j) {

        cout << "-------+";
    }
    cout << endl;

    cout << "min of column is: " << min_x << endl;
    cout << "max of row is:     " << max_y << endl;

    cout << "saddle point is ";
    if(isSaddlePoint()){
        cout << "exist.\n";
    }
    else{
        cout << "not exist.\n";
    }

    cout << "================================\n";
}

const vector<vector<double>> &Strategy::getMatrixForA() const {
    return Matrix;
}

// pivot operation
vector<vector<double>>& Strategy::createMatrixForB() {

    vector<vector<double>> result(Matrix.at(0).size(),vector<double>(Matrix.size()));
    for (size_t i = 0; i < Matrix.size(); ++i)
        for (size_t j = 0; j < Matrix.at(0).size(); ++j)
            result[j][i] = Matrix[i][j];

    MatrixB = result;
}

const vector<vector<double>> &Strategy::getMatrixForB() const {
    return MatrixB;
}

double Strategy::getMinX() const {
    return min_x;
}

double Strategy::getMaxY() const {
    return max_y;
}


class Simplex {

    vector<vector<double>> A;
    vector<double> B;
    vector<double> C;
    int size_x, size_y;
    double g;
    char gambler;


    vector<string> VARIABLES_Y;
    vector<string> VARIABLES_X;

    vector<vector<double>> full_matrix;


public:
    void print() {

        for (int i = 0; i < size_y; ++i) {
            cout << "+";
            for (int j = 0; j < size_x + 1; ++j) {
                cout << "-------+";
            }
            cout << endl;
            for (int j = 0; j < size_x; ++j) {
                std::cout << "|  " << roundf(A.at(i).at(j)*100)/100 << "\t";
            }
            cout << "|  " << roundf(B.at(i)*100)/100 << "\t|\n";
        }
        cout << "+";
        for (int j = 0; j < size_x + 1; ++j) {

            cout << "-------+";
        }
        cout << endl;
        for (int k = 0; k < size_x; ++k) {
            std::cout << "|  " << roundf(C.at(k)*100)/100 << "\t";
        }
        cout << "|\n";

        cout << "+";
        for (int j = 0; j < size_x; ++j) {

            cout << "-------+";
        }
        cout << endl;
    }
    void printFullMatrix(int y_size, int x_size) {

        cout << "FULL MATRIX IS:\n";

        for (int i = 0; i < y_size; ++i) {
            cout << "+";
            for (int j = 0; j < x_size; ++j) {
                cout << "-------+";
            }
            cout << endl;
            for (int j = 0; j < x_size; ++j) {
                std::cout << "|  " << roundf(full_matrix.at(i).at(j)*100)/100 << "\t";
            }
            cout << "|\n";
        }

        cout << "+";
        for (int j = 0; j < x_size; ++j) {
            cout << "-------+";
        }
        cout << endl;
    }
    void Solve() {
        if(gambler == 'A') {
            first_step_solve();
        }

        for (int i = 0; i < size_y; ++i) {
            printMatrixWithVariables();
            second_step_solve();
        }

        third_step_solve();

    }
    void first_step_solve() {

        for (int j = 0; j < size_y; ++j) {

            cout << "divide string at value in basics X" <<  j+5 << endl;

            double string_value_divider = full_matrix.at(j).at(j+size_x);
            for (int i = 0; i < full_matrix.at(0).size(); ++i) {
                if(full_matrix.at(j).at(i) != 0) {
                    full_matrix.at(j).at(i) = full_matrix.at(j).at(i) / string_value_divider;
                }
            }
            printFullMatrix(size_y,size_x+size_y+1);
        }
    }
    bool second_step_solve() {

        int column_index = find_max_or_min_index(); // find max or min from L(C) vector

        //get true string (divide B by a at index and then compare to choose min)
        int row_index = get_MIN_INDEX_BY_DIVIDING_AT_B(column_index);

        change_variables_index(column_index,row_index);


        divide_row_by_main_element(column_index,row_index);

        minus_main_row_from_all_other_rows(column_index,row_index);

        for (int i = 0; i < size_x; ++i) {
            C.at(i) = full_matrix.at(size_y).at(i);
        }
        cout << endl;

    }

    void printMatrix() {
        if(gambler == 'A') {
            printFullMatrix(size_y, size_x + size_y + 1);
        }
        else if (gambler == 'B') {
            printFullMatrix(size_y, size_x + size_y);
        }
    }


    void printMatrixWithVariables() {
        if(gambler == 'A') {
            printFullMatrixWithVariables(size_y +1, size_x + size_y + 1);
        }
        else if (gambler == 'B') {
            printFullMatrixWithVariables(size_y+1, size_x + size_y+1);
        }
    }


    void printFullMatrixWithVariables(int y_size, int x_size) {
        cout << "FULL MATRIX WITH VARIABLES IS:\n";

        cout << "+-------+";
        for (int j = 0; j < x_size ; ++j) {
            cout << "-----------+";
        }
        cout << endl;
        cout << "|       ";

        for (int d = 0; d < x_size-1; ++d){
            cout << "|  " <<  VARIABLES_X.at(d) << "\t\t";
        }
        cout << "|  " <<  "B" << "\t\t";
        cout << "|\n";

        for (int i = 0; i < y_size; ++i) {

            cout << "+-------+";
            for (int j = 0; j < x_size; ++j) {
                cout << "-----------+";
            }
            cout << endl;

            if(VARIABLES_Y.size() > i) {
                std::cout << "|  " << VARIABLES_Y.at(i) << 	"\t";
            } else {
                std::cout << "|  " << "L" << "\t";
            }


            for (int j = 0; j < x_size; ++j) {



                ostringstream strs;
                double a = roundf(full_matrix.at(i).at(j) * 100) / 100;
                strs << a;
                string str = strs.str();

                if (str.length() > 4) {
                    std::cout << "|  " << str << "\t";
                } else {
                    std::cout << "|  " << str << "\t\t";
                }

            }


            cout << "|\n";
        }

        cout << "+-------+";
        for (int j = 0; j < x_size; ++j) {
            cout << "-----------+";
        }
        cout << endl;
    }

    int find_max_or_min_index() {

        if(gambler == 'A'){
            //find max from L(C) vector
            int maxElementIndex = max_element(C.begin(),C.end()) - C.begin();
            return maxElementIndex;

        } else if (gambler == 'B') {
            //find min from L(C) vector
            int minElementIndex = min_element(C.begin(),C.end()) - C.begin();
            return minElementIndex;
        }
    }


    int get_MIN_INDEX_BY_DIVIDING_AT_B(int column) {
        vector<double> values;
        for (int i = 0; i < size_y; ++i) {
            values.emplace_back(abs(full_matrix.at(i).at(size_y+size_x)/full_matrix.at(i).at(column)));
        }

        cout << endl;

        int minElementIndex = std::min_element(values.begin(),values.end()) - values.begin();
        return minElementIndex;
    }


    void change_variables_index(int column, int row) {
        VARIABLES_Y.at(row) = VARIABLES_X.at(column);
    }


    void divide_row_by_main_element(int column, int row) {

        double element = full_matrix.at(row).at(column);
        for (int i = 0; i < size_x+size_y+1; ++i) {
            full_matrix.at(row).at(i) = full_matrix.at(row).at(i) / element;
        }
    }


    void minus_main_row_from_all_other_rows(int main_column, int main_row) {

        for (int i = 0; i < size_y+1; ++i) {
            if( i == main_row){
                continue;
            } else {

                double x_times = full_matrix.at(i).at(main_column) / full_matrix.at(main_row).at(main_column);

                for (int j = 0; j < size_x + size_y + 1; ++j) {
                    full_matrix.at(i).at(j) = full_matrix.at(i).at(j) - full_matrix.at(main_row).at(j) * x_times;
                }
            }
        }
    }

    void third_step_solve() {
        double sum = 0;
        for (int i = 0; i < VARIABLES_Y.size(); ++i) {
            int n = atoi(VARIABLES_Y.at(i).substr(1, 1).c_str());
            if(n < 5){
                sum += full_matrix.at(i).at(size_y+size_x);
                cout << "X" <<n <<": "<< full_matrix.at(i).at(size_y+size_x) << endl;
            }
        }
        if(gambler == 'A'){
            cout << "max win is " << 1/sum << endl;
        }
        else if (gambler == 'B') {
            cout << "min lose is " << 1/sum << endl;
        }
    }

    Simplex(vector<vector<double>> M, bool is_min, double min_win, double max_lose) {


        // FOR A = search min, FOR B = search max
        if (is_min) {
            gambler = 'A';
        } else {
            gambler = 'B';
        }

        A = M;
        size_y = M.size();
        size_x = M.at(0).size();

        for (int m = 0; m < size_y; ++m) {
            VARIABLES_Y.emplace_back("X" + to_string(m+1+size_x));
        }

        for (int l = 0; l < size_y+size_x; ++l) {
            VARIABLES_X.emplace_back("X" + to_string(l+1));
        }

        int c = 0;

        if (gambler == 'A') {
            c = -1;
        } else if (gambler == 'B') {
            c = 1;
        }


        // <= или  >= система

        for (int i = 0; i < size_y; ++i) {
            B.emplace_back(c);
        }


        for (int j = 0; j < size_x; ++j) {
            C.emplace_back(-1);
        }

        if (gambler == 'A') {
            g = (1 / min_win);
        } else {
            g = (1 / max_lose);
        }

        if (gambler == 'A') {
            vector<vector<double>> v(size_y + 1, vector<double>(size_x + size_y + 1, 0));
            full_matrix = v;
        } else
        if (gambler == 'B') {
            vector<vector<double>> v(size_y + 1, vector<double>(size_x + size_y+1, 0));
            full_matrix = v;
        }



        {
            for (int k = 0; k < size_y; ++k) {
                for (int i = 0; i < size_x; ++i) {
                    full_matrix.at(k).at(i) = M.at(k).at(i);
                }
            }


            for (int j = 0; j < size_y; ++j) {
                full_matrix.at(j).at(j+size_x) = B.at(j);
                if(gambler == 'A') {
                    full_matrix.at(j).at(size_x + size_y) = -1;
                }
            }

            for (int l = 0; l < size_x; ++l) {
                full_matrix.at(size_y).at(l) = C.at(l);
            }

            for(int i = 0; i < size_y; ++i) {
                if(gambler == 'A') {
                    full_matrix.at(i).at(size_x + size_y) = -B.at(i);
                } else if (gambler == 'B') {
                    full_matrix.at(i).at(size_x + size_y) = B.at(i);
                }
            }
        }
    }
};


int main() {

    vector<vector<double>> strategy_a1;
    strategy_a1.push_back({{ 6, 18,  6, 15}});
    strategy_a1.push_back({{17,  8, 13, 14}});
    strategy_a1.push_back({{16, 16, 18, 2}});
    strategy_a1.push_back({{18,  8,  4, 18}});
    strategy_a1.push_back({{15,  8,  3, 19}});

    Strategy A(strategy_a1);
    A.print();
    auto M = A.getMatrixForA();
    Simplex K(M, true,A.getMaxY(),A.getMinX());
    K.print();
    K.printMatrix();
    K.printMatrixWithVariables();
    K.Solve();
    auto M1 = A.getMatrixForB();
    Simplex K1(M1, false, A.getMaxY(), A.getMinX());
    K1.print();
    K1.printMatrix();
    K1.printMatrixWithVariables();
    K1.Solve();

    return 0;
}
