# matrix-3-all-21

## unsat ==> sat (watched inc)

## clauses:

```z3-trace
0
true

1
!(x20 < 0)

2
!(x41 < 0)

3
!(x4 < 0)

4
!(x21 < 0)

5
!(x2 < 0)

6
!(x36 < 0)

7
2 x2 + x4 > 0

8
!(2 x2 + x4 < 0)

9
!(x2 + 2 < 0)

10
x21 x2 - 2 x2 - x4 + x20 < 0

11
!(x21 x2 - 2 x2 - x4 + x20 > 0)

12
x21 x36 - 2 x36 - x4 + x20 < 0

13
!(x21 x36 - 2 x36 - x4 + x20 > 0)

14
!(x21 - 2 > 0)

15
!(- x21 + x41 x4 - x20 x41 + 2 < 0)

16
x21 x36 + x20 > 0

17
!(x21 x36 + x20 < 0)

18
!(x21 + x20 x41 < 0)

19
x4 - x20 < 0

20
!(x4 - x20 > 0)

21
!(x21 - 2 < 0)
```

## watched clauses:

```z3-trace
var 0 x20
clause 18
!(x21 + x20 x41 < 0)
clause 19
x4 - x20 < 0
clause 20
!(x4 - x20 > 0)

var 1 x41
clause 15
!(- x21 + x41 x4 - x20 x41 + 2 < 0)

var 2 x4
clause 7
2 x2 + x4 > 0
clause 8
!(2 x2 + x4 < 0)
clause 19
x4 - x20 < 0
clause 20
!(x4 - x20 > 0)

var 3 x21
clause 10
x21 x2 - 2 x2 - x4 + x20 < 0
clause 11
!(x21 x2 - 2 x2 - x4 + x20 > 0)
clause 12
x21 x36 - 2 x36 - x4 + x20 < 0
clause 13
!(x21 x36 - 2 x36 - x4 + x20 > 0)
clause 15
!(- x21 + x41 x4 - x20 x41 + 2 < 0)
clause 16
x21 x36 + x20 > 0
clause 17
!(x21 x36 + x20 < 0)
clause 18
!(x21 + x20 x41 < 0)

var 4 x2
clause 7
2 x2 + x4 > 0
clause 8
!(2 x2 + x4 < 0)
clause 10
x21 x2 - 2 x2 - x4 + x20 < 0
clause 11
!(x21 x2 - 2 x2 - x4 + x20 > 0)

var 5 x36
clause 12
x21 x36 - 2 x36 - x4 + x20 < 0
clause 13
!(x21 x36 - 2 x36 - x4 + x20 > 0)
clause 16
x21 x36 + x20 > 0
clause 17
!(x21 x36 + x20 < 0)
```















