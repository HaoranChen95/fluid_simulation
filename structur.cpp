// struct.cpp -- a simple structure

# include <iostream>
# include <string.h>

using namespace std;

struct inflatable // 定义结构，用于同一变量存储不同值，单个值可以改变
{
    string name;
    float volume;
    double price;
};

union one4all //共用体，可以存储不同类型数据可是，只能存储一个值，但并不智能
{
    int int_val;
    long long_val;
    double double_val;
};

enum specrom {r=1,g=2,b=3};

int main()
{
    //使用结构
    struct inflatable guest =
    {
        "Glorious Gloria",
        1.88,
        29.99
    };

    struct inflatable pal =
    {
        "Audacious Arthur",
        3.12,
        32.99
    };
    cout << guest.name << endl;

    guest.volume = 100;
    guest.name = "new name";

    cout << guest.volume << endl;

    //使用共用体
    one4all numb;

    numb.double_val = 1.0;
    cout << numb.double_val <<endl;

    //使用枚举
    
    int color = r; //去调用枚举中的值
    cout << color << endl; //取值范围是2的n次幂-1
    
    return 0;
}
