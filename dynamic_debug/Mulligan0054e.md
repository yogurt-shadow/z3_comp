# MulliganEconomicsModel0054e

## Sat ==> Unsat (random)

```z3-trace
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b1 x0 < 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b2 (x2)(x1) > 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b3 (x1)(- x3 + x1 x2) > 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b4 x3 = 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b5 x1 > 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b6 x3 > 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b7 x4 - 1 < 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b8 x4 < 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b9 x5 + 1 > 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b10 x0 x1^3 x6 - x3^2 - x1^3 = 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b11 x1 x7 + x1 x6 + x3 - x1 x2 + x1 = 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b12 x5 + 1 = 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b13 - x7 + x5 x6 + x6 = 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b14 x8 = 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b15 x6 + 1 > 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b16 x6 < 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b17 x7 < 0
------------------------------------------------
-------- [nlsat_verbose] mk_ineq_atom ../src/nlsat/nlsat_solver.cpp:1087 ---------
create: b18 (x1)(- x3 + x1 x2 - x1) < 0
```

```z3-trace
v10 < 0
(v8)(v2) > 0
(v2)(- v1 + v2 v8) > 0
!(v1 = 0)
v2 > 0
v1 > 0
v4 - 1 < 0
!(v4 < 0)
v3 + 1 > 0
v10 < 0
v10 v2^3 v5 - v1^2 - v2^3 = 0
v2 v7 + v2 v5 + v1 - v2 v8 + v2 = 0
v3 + 1 = 0 or - v7 + v3 v5 + v5 = 0
v5 + 1 > 0
v5 < 0
v7 < 0
(v2)(- v1 + v2 v8 - v2) < 0
```



| var  | index | assignment |
| ---- | ----- | ---------- |
| v1   | 3     | -1         |
| v5   | 6     | -0.5       |
| v8   | 2     | 0.125      |
| v6   | 8     | 0.5        |
| v7   | 7     | -1         |
| v10  | 0     | -1         |
| v4   | 4     | 0.5        |
| v3   | 5     | 1          |
| v2   | 1     | conflict   |

+ conflict:

$$
v_8 v_2>0\\
v_{10} v_2^3v_5-v_1^2-v_2^3=0
$$

+ new valid clause:

$$
!(v_8 > 0)||!(v_{10}>root[1](v_{10} v_5-1))||!(v_5<0)|| !(v_1>0)
$$

+ backtrack to v_10
| var  | index | assignment |
| ---- | ----- | ---------- |
| v1   | 3     | -1         |
| v5   | 6     | -0.5       |
| v8   | 2     | 0.125      |
| v6   | 8     | 0.5        |
| v7   | 7     | -1         |
| v10  | 0     | -1 ==> -3  |
| v2   | 1     | conflict   |

+ conflict:

$$
v_{10} v_2^3v_5-v_1^2-v_2^3=0\\
v_2(-v_1+v_2v_8)>0
$$

+ new valid clause:

$$
!(v_{10} < root[1](v_{10} v_1 v_5 - v_1 - v_8^3)) || !(v_8 > 0) \\
|| !(v_5 < 0) || !(v_1 > 0)
$$

+ backtrack to v_10
| var  | index | assignment          |
| ---- | ----- | ------------------- |
| v1   | 3     | -1                  |
| v5   | 6     | -0.5                |
| v8   | 2     | 0.125               |
| v6   | 8     | 0.5                 |
| v7   | 7     | -1                  |
| v10  | 0     | -3 ==> -2.001953125 |
| v4   | 4     | 0.5                 |
| v3   | 5     | 1                   |
| v2   | 1     | conflict            |

+ conflict:

$$
v_2 v_7 + v_2 v_5 + v_1 - v_2 v_8 + v_2 = 0\\
v_{10}v_2^3-v_1^2-v_2^3=0
$$

+ new valid clause:

$$
v1^2 v7^3 + 3 v1^2 v5 v7^2 - 3 v8 v1^2 v7^2 + 3 v1^2 v7^2 + 3 v1^2 v5^2 v7 - 6 v8 v1^2 v5 v7 + 6 v1^2 v5 v7 + 3 v8^2 v1^2 v7 - 6 v8 v1^2 v7 + 3 v1^2 v7 + v1^2 v5^3 - 3 v8 v1^2 v5^2 + 3 v1^2 v5^2 + v10 v1^3 v5 + 3 v8^2 v1^2 v5 - 6 v8 v1^2 v5 + 3 v1^2 v5 - v1^3 - v8^3 v1^2 + 3 v8^2 v1^2 - 3 v8 v1^2 + v1^2 = 0 \\
|| v7 + v5 - v8 + 1 = 0  \\
|| 
!(v7 + v5 - v8 + 1 < 0)
$$

+ backtrack to v_10: conflict



















