'''
Author: wzh
LastEditors: wzh
Date: 2022-05-06 10:29:12
Description: file content
FilePath: \dynamic_debug\test.py
LastEditTime: 2022-05-06 10:48:10
'''
from fractions import Fraction

if __name__ == "__main__":
    v1 = Fraction(1, 4)
    v8 = Fraction(1, 2)
    v10 = -2
    v5 = -Fraction(17, 32)
    v7 = -Fraction(7, 32)
    res = v1*v1 *v7*v7*v7 + 3 *v1*v1* v5 *v7*v7 - 3 *v8 *v1*v1 *v7*v7 + 3 *v1*v1* v7*v7 + \
3 *v1*v1* v5*v5* v7 - 6 *v8 *v1*v1* v5* v7 + 6 *v1*v1* v5 *v7 + 3* v8*v8* v1*v1* v7 - 6 *v8 *v1*v1* v7 + \
3 *v1*v1* v7 + v1*v1* v5*v5*v5 - 3 *v8 *v1*v1* v5*v5 + 3* v1*v1* v5*v5 + v10 *v1*v1*v1* v5 + \
3* v8*v8* v1*v1* v5 - 6 *v8* v1*v1* v5 + 3 *v1*v1* v5 - v1*v1*v1 - v8*v8*v8* v1*v1 + 3 *v8*v8* v1*v1 \
- 3 *v8* v1*v1 + v1*v1
    print(res)
    v2 = 1
    res2 = v2 *v7 + v2* v5 + v1 - v2 *v8 + v2
    print(res2)

