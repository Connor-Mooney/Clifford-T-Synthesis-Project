from Operator import Operator
# THIS FILE WILL HAVE THE MAIN MACHINERY WORKING
Patterns = [[],[],[],[],[]]
alls = []
for i in range(6):
    for j in range(i+1,6):
        if (i,j) in ((0,1),(1,2),(3,4),(4,5),(1,0),(2,1),(4,3),(5,4)):
            sign = 1
        else:
            sign = -1
        name = f'T[{i+1},{j+1}]'
        k = 1
        mat = [np.zeros((6,6)),np.identity(6)]
        mat[0][i,i] = 1
        mat[0][i,j] = -sign
        mat[0][j,i] = sign
        mat[0][j,j] = 1
        mat[1][i,i] = 0
        mat[1][j,j] = 0
        op = Operator(k,mat,name)
        Patterns[0].append(op)
        alls.append(op.permclass)
print(len(Patterns[0]))

for i in range(1,5):
  for t in Patterns[i-1]:
    for op in Patterns[i-1]:
     prod = t.dot(op)
     if prod.permclass not in alls:
        Patterns[i-1].append(prod)
        alls.append(prod.permclass)

