#include <iostream>
using namespace std;
int location(double ***array,int row, int col, int depth){
    int depth_num = sizeof(array[0][0]);
    int row_num = sizeof(array[0])/(depth_num/sizeof(double));
    int col_num = sizeof(array)/(depth_num/sizeof(double))*(row_num/sizeof(double));
    int locale = 2000 + row * col_num * depth_num + col * depth_num + depth;
    return locale;
}
double *** array_constructor(int rows , int cols, int depth){
    double ***array;
    array = new double**[rows];
    for(int i = 0; i < 5; i++){
        array[i] = new double*[cols];
        for(int j = 0; j < 5; j++){
            array[i][j] = new double[depth];
        }
    }
    return array;
}
int main(){
    double ***array = array_constructor(3,2,5);
    cout << location(array,1,1,4);
    cout << '\n';
    return 0;
}