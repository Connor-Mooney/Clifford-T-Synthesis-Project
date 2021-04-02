class Operator:
  """This class represents a matrix in SO6(Z[1/sqrt(2)]), storing the name of the operator, i.e. a string of T operators, its actual matrix, its residue, 
  and its column-permutation-invariant class"""
  def __init__(self, mat,name):
    self.matrix = mat
    self.name = name
    self.residue = self.gen_res()
    self.permclass = self.column_perm_class()
  
  def gen_res(self):
    """Generates the residue of this operator"""
    return None
  def column_perm_class(self):
    """Generates the permutation class of this operator"""
    return None
  def dot(self, op2):
    """multiplies this operator with another"""
    return None
    
