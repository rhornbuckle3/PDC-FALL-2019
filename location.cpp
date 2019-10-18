//Russell Hornbuckle
//2019
//location of an element in memory in a 3 dimensional array that is hypothetically initialized at 2000
#include <iostream>
using namespace std;
int location(double ***array,int row, int col, int depth, int row_num, int col_num,int depth_num){
    //can't calculate dimensions with a pointer array unfortunately
    int locale = 2000 + (row * col_num * depth_num + col * depth_num + depth)*sizeof(array);
    return locale;
}
double *** array_constructor(int rows , int cols, int depth){
    double ***array;
    array = new double**[rows];
    for(int i = 0; i < rows; i++){
        array[i] = new double*[cols];
        for(int j = 0; j < cols; j++){
            array[i][j] = new double[depth];
            //could fill it here with another loop
        }
    }
    return array;
}
int main(){
    double ***array = array_constructor(3,2,5);
    cout << location(array,1,1,4,3,2,5);
    cout << '\n';
    return 0;
}