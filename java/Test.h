/*************************************************************************
	> File Name: Test.h
	> Author: ma6174
	> Mail: ma6174@163.com 
	> Created Time: Sat Jul 20 13:01:15 2019
 ************************************************************************/

#include<iostream>

using namespace std;


class Test{
	public:
	int val = 0;
	Test(){
		std::cout<<"test constructor"<<std::endl;
	}
	~Test(){
		std::cout<<"test deconstructor"<<std::endl;
	}

};
