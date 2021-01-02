# include <iostream>
//地址运算符&, 指针运算符运算符*

using namespace std;

int main()
{
    int int_val = 5;
    int * pointer_2_int; //定义指针，多变量定义时，要对每一个变量加上*，系统只分配了存储地址的内存，还没分配存储值的内存
    pointer_2_int = &int_val; //给指针赋值 指针需要初始化

    cout << &int_val << endl; //获得地址
    cout << pointer_2_int << endl; //指针值
    cout << int_val <<endl; //参数值
    cout << *pointer_2_int << endl;//指针对应的参数值

    //使用new动态分配内存
    int *pn = new int;
    *pn = 1000;
    cout << *pn << "\n" << pn <<endl;
    delete pn; //释放内存

    //使用new创建动态数组
    int *psome = new int [10]; //psome是第一个值的地址, 需要确定值的范围
    psome[0] =0;
    psome[1] =1;
    psome[2] =2;
    cout << psome[1] << endl;
    psome++;
    cout << psome[1] << endl;
    psome--;
    delete [] psome; 

    return 0;
}