"""building WAPB functions from https://link.springer.com/article/10.1007/s10623-018-0579-x """

def oj(n,j):
  """Return the cardinality of the cyclotomic class of 2 modulo 2 n − 1 containing j"""
  return len(Zmod(2^n-1).cyclotomic_cosets(2, cosets=[j])[0])
  
def my_tr_1(k,x,F):
  """Returns the sum_{i=0}^{k-1} x^(2^i) over the field F"""
  return F(sum([x^(2^i) for i in range(k)]))
  
def random_LM(n):
  """Returns a function in n variables sampled uniformly at random from LM family

  Args:
      n (int): number of variables

  Returns:
      BooleanFunction: LM-Boolean function in n variables 
  """

  F=GF(2^n, 'z')
  z=F.gen()
  L0=[l[0] for l in Zmod(2^n-1).cyclotomic_cosets(2) if l[0]!=0]
  FR.<x>=PolynomialRing(F) 
  SF=F.subfields()
  F2,m2=SF[1]
  b=F2.gen()
  L=[my_tr_1(oj(n,j),m2(b^(ZZ.random_element(1,3)))*x^j,FR) for j in L0 ] 
  fp=sum(L)
  T=[fp(a) for a in F]
  f=BooleanFunction(T)
  return f
  

def all_LM(n):
  """Returns a function in n variables sampled uniformly at random from LM family
  
  Args:
      n (int): number of variables

  Returns:
      list[BooleanFunction]:  LM family in n variables 
  
  Example:

      LM4=all_LM(4)
      for f in LM4:
        g=construction_2(4,f)
        table_fun(g)

  """
  F=GF(2^n, 'z')
  z=F.gen()
  Gamma0=[l[0] for l in Zmod(2^n-1).cyclotomic_cosets(2) if l[0]!=0]
  FR.<x>=PolynomialRing(F) 
  SF=F.subfields()
  F2,m2=SF[1]
  b=F2.gen()
  gamma0=len(Gamma0)
  #print gamma0
  alln=[]
  for v in prodit(range(1,3), repeat=gamma0):
    fp=FR(0)
    for i in range(gamma0):
      j=Gamma0[i]
      fp+=my_tr_1(oj(n,j),m2(b^(v[i]))*x^j,FR)
    T=[fp(a) for a in F]
    f=BooleanFunction(T)
    #print is_WPB(f),
    #print [NLk(k,f) for k in [2,3,4]]
    alln+=[f]
  return alln
  

def LM8():
  """Returns a function of LM function with nlk list [6,21,27,22,9] in 8 variables 
     This is l in Table 4 of https://eprint.iacr.org/2022/1434.pdf
  """  
  n=8
  F=GF(2^n, 'z')
  z=F.gen()
  L0=[l[0] for l in Zmod(2^n-1).cyclotomic_cosets(2) if l[0]!=0]
  v=ones_matrix(len(L0),1)
  v[-4]=2
  FR.<x>=PolynomialRing(F) 
  SF=F.subfields()
  F2,m2=SF[1]
  b=F2.gen()
  fp=FR(0)
  for i in range(len(L0)):
    j=L0[i]
    fp+=my_tr_1(oj(n,j),m2(b^(v[i][0]))*x^j,FR)
  T=[fp(a) for a in F]
  f=BooleanFunction(T)
  return f 

def LM8_():
  """Returns a WPB function, derived from an LM function, with nlk list [9,22,27,22,9] in 8 variables 
     This is l0 in Table 4 of https://eprint.iacr.org/2022/1434.pdf
  """    
  f=LM8()
  n=8
  nf= rev_fun(f)
  L=[NLk(k,f) for k in range(0,9)]
  Ln=[NLk(k,nf) for k in range(0,9)]
  h=Ikn(0,n)*f
  for k in range(1,9):
    if L[k]>=Ln[k]:
      h+=Ikn(k,n)*f
    else:
      h+=Ikn(k,n)*nf    
  return h