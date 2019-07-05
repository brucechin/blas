#include<iostream>
#include<sys/time.h>
#include<cblas.h>
#include<vector>
#include<unistd.h>
#include<cstdlib>
#include<string>


class LogicMatrix{
private:

    int nrow;
    int ncol;
    bool* value;

public:

    LogicMatrix(){

    }

    LogicMatrix(int n, int m){

    }

    LogicMatrix(bool* mat){

    }

    LogicMatrix(LogicMatrix logic){

    }

    int getNRow(){
        return this.nrow;
    }

    int getNCol(){
        return this.ncol;
    }

    bool* getLogicMatrix(){
        return this.value
    }

    bool getElement(int i, int j){
        return this.value[i][j];
    }

    bool getRowVector(int i){

    }

    double getRowVectorAsDouble(int i){

    }

    void setElement(int i, int j, bool x){

    }

    //set all value to x
    void setValue(bool x){

    }

    void copyTo(LogicMatrix mat){

    }




}