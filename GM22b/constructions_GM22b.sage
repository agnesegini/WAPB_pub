"""
  Constructions described in [GM22b] https://eprint.iacr.org/2022/1434
"""


def Sn(n):
  """Sn={k ∈ [1, n/2[ : {n choose k-1} = {n choose k} = 1 mod 2} See Remark 1 [GM22b]"""  
  if (n%2 == 0) :  return []
  return list(filter(lambda e : (((binomial(n,e-1)%2)==1) and (binomial(n,e)%2)==1) , range(1,ceil(n/2))))

def construction_1_GM22b(n,f,g):
  """Secondary construction from two n-variables SWAPB functions to one n+1 SWAPB function.
     This is Construction 1 from [GM22b].

  Args:
      n (int): number of variables input functions
      f (BooleanFunction): n-variable SWAPB function
      g (BooleanFunction): n-variable SWAPB functions

  Returns:
      BooleanFunction: n+1-variables SWAPB function
  """  
  for k in Sn(n):
    g=g+Ikn(k-1,n)+Ikn(n-k,n)
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n+1)])
  ga=R(g.algebraic_normal_form())
  fa=R(f.algebraic_normal_form())
  newf=R(fa)*R('1 + x%d' %n)  +R(ga)*R('x%d' %n)
  print(newf)
  return BooleanFunction(newf)


def construction_2_GM22b(t,f, verbose=False):
  """Secondary construction from a n-variables SWAPB function to one n+t SWAPB function.
  Let n, t ∈ N ∗ f a n-variable SWAPB functions.
Output: g an 

  Args:
      t (int): number of new variables
      f (int): n-variables SWAPB function
      verbose (bool, optional): printing details. Defaults to False.

  Returns:
      BooleanFunction: (n + t)-variable SWAPB function.
  """  
  n=f.nvariables()
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n+t)])
  g=f
  for i in range(1,t+1):
    if verbose: print(i)
    u=n+i-1
    S=Sn(u)
    h=BooleanFunction(u)
    for k in S:
      h+=Ikn(k-1,u)+Ikn(u-k,u)
      
    g=BooleanFunction(R(g.algebraic_normal_form()))+BooleanFunction(R('x%d' %u)*R(h.algebraic_normal_form()))
  return g


def no_int_bin_n(d,n):
  """Retuns the binary decomposition of the integer d xor 1 """
  d=bin(d)[2:]
  ld=len(d)
  assert (n-ld)>=0, ld
  d=(n-ld)*'0'+d
  return vector(ZZ,vector(GF(2),d)+vector(ones_matrix(GF(2),n,1)))


def no_int_n(d,n):
  """Retuns the integer d xor 1"""
  s=no_int_bin_n(d,n)
  return sum([s[i]*2^(n-i-1) for i in range(n)])

def rev_fun(f):
  """Retuns f(x xor 1)"""
  n=f.nvariables()
  T=[f(no_int_n(d,n)) for d in range(2^n)]
  return BooleanFunction(T)   
  
