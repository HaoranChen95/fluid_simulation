// string.app -- storing strings in an array

# include <iostream>
# include <cstring> // for the strlen() function

int main()
{
    using namespace std;
    const int Size =15;
    char name1[Size] ;
    char name2[Size] = "C++owboy";
    cout << "Howdy! I'm " << name2;
    cout << "! What's your name?\n";
    //cin.getline(name1, Size);  cin 读取的范围到空格会停止，如果要获得多个单词应当使用cin.getline
    // 或者用
    cin.get(name1, Size).get();
    cout << "Well, " << name1 << ", your name has ";
    cout << strlen(name1) << " letters and is stored\n";
    cout << "in an array of " << sizeof(name1) << " bytes\n";
    cout << "Your initial is " << name1[0] << endl;
    name2[3] = '\0' ;
    cout << "Here the first 3 characters of my name: ";
    cout << name2 <<endl;
    return 0;
}