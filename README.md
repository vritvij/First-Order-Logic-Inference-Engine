# First-Order-Logic-Inference-Engine
A C++ program that can tell statements to and ask a query from a first order logic based knowledge base. Inserting statements into the knowledge base required converting them to CNF form before transforming them to the knowledge base's internal data structure which is optimized for their retrieval and processing. The resolution inference algorithm is used when the knowledge base is asked a query.

##Input Format: Input.txt
```
NQ: (int) The number of queries that you want to ask the knowledge base
Queries: (string[NQ]) The queries
NS: (int) The number of first order logic statements in the knowledge base
Statements: (string[NS]) The first order logic statements
```

##Output Format: Output.txt
```
Outputs: (string[NQ]) The output of the queries on the knowledge base (TRUE/FALSE)
```

##Examples
```
2
Ancestor(Liz,Billy)
Ancestor(Liz,Bob)
6
Mother(Liz,Charley)
Father(Charley,Billy)
((~Mother(x,y)) | Parent(x,y))
((~Father(x,y)) | Parent(x,y))
((~Parent(x,y)) | Ancestor(x,y))
((~(Parent(x,y) & Ancestor(y,z))) | Ancestor(x,z))
```
```
6
F(Bob)
H(John)
~H(Alice)
~H(John)
G(Bob)
G(Tom)
14
(A(x) => H(x))
(D(x,y) => (~H(y)))
((B(x,y) & C(x,y)) => A(x))
B(John,Alice)
B(John,Bob)
((D(x,y) & Q(y)) => C(x,y))
D(John,Alice)
Q(Bob)
D(John,Bob)
(F(x) => G(x))
(G(x) => H(x))
(H(x) => F(x))
(R(x) => H(x))
R(Tom)
```
```
2
Student(Alice)
Student(Bob)
1
(Student(a) | (Student(b) | (Student(c) | (Student(d) | (Student(e) | Student(f))))))
```
