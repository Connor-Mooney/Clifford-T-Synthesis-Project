import numpy as np

class Operator:
  """This class represents a matrix in SO6(Z[1/sqrt(2)]), storing the name of the operator, i.e. a string of T operators, its actual matrix, its residue, 
  and its column-permutation-invariant class"""
  def __init__(self, LDE, mat,name):
    self.matrix = mat
    self.LDE = LDE
    self.name = name
    self.residue = self.gen_res()
    self.permclass = self.column_perm_class()
  
  def gen_res(self):
    """Generates the residue of this operator"""
    return [self.matrix[0] % 2, self.matrix[1] % 2]
  
  def transpose(mat):
    """returns the transpose of a matrix"""
    vOut = []
    for i in range(len(mat[0])):
        vOut.append([mat[j][i] for j in range(len(mat))])
    return vOut

  
  def column_perm_class(self):
    """Generates the permutation class of this operator"""
    permclass = [[None, None, None, None, None, None], [None, None, None, None, None, None], \
                 [None, None, None, None, None, None], [None, None, None, None, None, None],\
                 [None, None, None, None, None, None], [None, None, None, None, None, None]]
    
    for i in range(self.matrix[0].length):
      for j in range(self.matrix[0][0].length):
        permclass[i][j] = (self.matrix[0][i,j], self.matrix[1][i,j])
   
    return [self.tcount, transpose(sorted(transpose(permclass)))]
    
 
  def dot(self, op2):
    """multiplies this operator with another"""
    LDENew = self.LDE + op2.LDE
    mat0 = np.matmul(self.matrix[0],op2.matrix[0]) + 2*np.matmul(self.matrix[1],op2.matrix[1])
    mat1 = np.matmul(self.matrix[0],op2.matrix[1]) + np.matmul(self.matrix[1],op2.matrix[0])
    newmat = [mat0, mat1]
    newname = self.name + " " + op2.name
    while np.all(newmat[0] % 2 == 0):
      # NEED TO REDUCE THIS DOWN TO MOST SIMPLIFIED MATRIX
      mat0 = np.copy(newmat[1])
      mat1 = np.copy(newmat[0])/2
      newmat = [mat0, mat1]
      LDENew = LDENew - 1
    return Operator(LDE, newmat, newname) 
    
