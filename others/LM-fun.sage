## building  CMR WAPB functions from https://link.springer.com/article/10.1007/s10623-018-0579-x

#return the cardinality of the cyclotomic class of 2 modulo 2 n − 1 containing j
def oj(n,j):
  return len(Zmod(2^n-1).cyclotomic_cosets(2, cosets=[j])[0])
  
#return the sum_{i=0}^{k-1} x^(2^i) over the field F
def my_tr_1(k,x,F):
  return F(sum([x^(2^i) for i in range(k)]))
  
#return random LM  WPB function in n variables 
def random_LM(n):
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
  
#return all the LM  WPB function in n variables  
def all_LM(n):
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
  
#return a LM function l of Table 4 with nlk list [6,21,27,22,9]  in 8 variables 
def LM8():
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
  
#return a LM function l0 of Table 4 with nlk list [9,22,27,22,9]  in 8 variables 
def LM8_():
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
  
"""
LM4=all_LM(4)
for f in LM4:
  g=construction_2(4,f)
  table_fun(g)
""" 
