/*
 * @Author: wzh
 * @LastEditors: wzh
 * @Date: 2022-05-06 10:16:52
 * @Description: file content
 * @FilePath: \Desktop\test.cpp
 * @LastEditTime: 2022-05-06 10:22:04
 */
#include <iostream>

using namespace std;

int main(){
double v1 = 0.25, v7 = -0.21875, v5 = -0.53215, v10 = -2, v8 = 0.5;
double res = v1*v1 *v7*v7*v7 + 3 *v1*v1* v5 *v7*v7 - 3 *v8 *v1*v1 *v7*v7 + 3 *v1*v1* v7*v7 + 
3 *v1*v1* v5*v5* v7 - 6 *v8 *v1*v1* v5* v7 + 6 *v1*v1* v5 *v7 + 3* v8*v8* v1*v1* v7 - 6 *v8 *v1*v1* v7 + 
3 *v1*v1* v7 + v1*v1* v5*v5*v5 - 3 *v8 *v1*v1* v5*v5 + 3* v1*v1* v5*v5 + v10 *v1*v1*v1* v5 + 
3* v8*v8* v1*v1* v5 - 6 *v8* v1*v1* v5 + 3 *v1*v1* v5 - v1*v1*v1 - v8*v8*v8* v1*v1 + 3 *v8*v8* v1*v1 
- 3 *v8* v1*v1 + v1*v1;
cout << res << endl;
    return 0;
}