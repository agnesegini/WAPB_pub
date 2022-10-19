
def Sn(n):
  return list(filter(lambda e : (((binomial(n,e-1)%2)==1) and (binomial(n,e)%2)==1) , range(1,ceil(n/2))))

def construction_1(n,f,g):
  for k in Sn(n):
    g=g+Ikn(k-1,n)+Ikn(n-k,n)
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n+1)])
  ga=R(g.algebraic_normal_form())
  fa=R(f.algebraic_normal_form())
  newf=R(fa)*R('1 + x%d' %n)  +R(ga)*R('x%d' %n)
  print(newf)
  return BooleanFunction(newf)


def construction_2(t,f, verbose=False):
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

#retuns the binary decomposition of the integer d xor 1 
def no_int_bin_n(d,n):
  d=bin(d)[2:]
  ld=len(d)
  assert (n-ld)>=0, ld
  d=(n-ld)*'0'+d
  return vector(ZZ,vector(GF(2),d)+vector(ones_matrix(GF(2),n,1)))


#retuns the integer d xor 1 
def no_int_n(d,n):
  s=no_int_bin_n(d,n)
  return sum([s[i]*2^(n-i-1) for i in range(n)])

#retuns f(x xor 1) 
def rev_fun(f):
  n=f.nvariables()
  T=[f(no_int_n(d,n)) for d in range(2^n)]
  return BooleanFunction(T)   
  
