import random

def elem_2(n):
  """Symmetric function of degree 2 in n-variables

  Args:
      n (int): number of variables

  Returns:
      BooleanPolynomial
  """  
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  f=R(0)
  for i in range(n):
    for j in range(i):
     f+=R('x%d*x%d' %(i,j))
  return f
      
def elem_symm(d,n):
  """Symmetric function of degree d in n-variables.

  Args:
      d (int): degree
      n (int): number of variables

  Returns:
      BooleanPolynomial
  """  
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  f=R(0)
  for ind in ithp(n,d):
     s='x%d'+(d-1)*'*x%d'
     f+=R(s %ind)
  return f      
      
def ell_half(n):
  """Sum of the first n/2 variables of BooleanPolynomialRing in n variables.

  Args:
      n (int): number of variables

  Returns:
      BooleanPolynomial
  """  
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  return sum([R('x%d' %i) for i in range(0,n//2)])    

def magic(n):
  """Boolean function give by symmetric function of degree 2 in n-variables and first n/2 variables.

  Args:
      n (int): number of variables

  Returns:
      BooleanFunction
  """
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  f=sum([R('x%d' %i) for i in range(0,n//2)])
  for i in range(n):
    for j in range(i):
     f+=R('x%d*x%d' %(i,j))
  return BooleanFunction(f)

def comp(f,n):
  """Produces the sum of the absolute values of the Walsh transforms at 0 restricted to all the slices
  """  
  return sum([abs(walsh_tr_0(f,Ekn(i,n))) for i in range(n+1)])

def NPB(f):
  """Returns the Non Perfect Balancedness of a boolean function (Definition 11 [GM22c])

  Args:
      f (BooleanFunction): Boolean function

  Returns:
      int : NPB(f)
  """   
  n=f.nvariables()
  o=(2-walsh_tr_0(f,Ekn(0,n))+walsh_tr_0(f,Ekn(n,n)))
  t=sum([abs(walsh_tr_0(f,Ekn(i,n))) for i in range(1,n)])
  return (t+o)/2


def count_GM22c(f):
   """Returns the number of possible distinct outputs of construction1_GM22c(f) """
   #f=magic(n)
   n=f.nvariables()
   L=[walsh_tr_0(f,Ekn(i,n)) for i in range(n+1)]
   print(L)
   t=1
   for k in range(1,n):
      t=t*binomial(binomial(n,k)/2+L[k]/2,L[k])
   print(t,log(t,2).n())
   return t

def construction1_GM22c(seed,verbose=False):
  """
  Construction 1 from [GM22c]  to obtain a WPB function from any 2^m -variable Boolean function seed.

  Args:
      seed (BooleanFunction): input function
      verbose (bool, optional): Prints info at every loop. Defaults to False.

  Returns:
      (BooleanFunction,dict): The first element is the the output function. The second one is a
      dictionary D such that D[k] is the list of points in the kth slice whose output bit was flipped.
  """
  n=seed.nvariables()
  FT=seed.truth_table('int')
  N=list(FT)
  N[0]=0
  N[-1]=1
  wpb, t =is_WPB(BooleanFunction(N))
  flipBooleanFunctionped={ 0: [0], n: [0] }
  while not wpb :
    E=Ekn(t,n)
    o=walsh_tr_0(seed,E)//2
    S=suppk(seed,t)
    if verbose: print(t,o)
    if o>0:
      CS = [x for x in E if x not in S] 
      O=sample(CS, o)
      for i in O:
        N[i]=1
    else:
      O=sample(S, -o)
      for i in O:
        N[i]=0
    flipped[t]=O  
    y=BooleanFunction(N)
    wpb, t =is_WPB(y) 
  nly=y.nonlinearity()
  return  y,flipped,nly

def upmagic(n):
  """Construction 1 from [GM22c] applied with magic(n) as input (bit-flipping deterministic)"""
  f=magic(n)
  FT=f.truth_table('int')
  N=list(FT)
  N[-1]=1
  wpb, t =is_WPB(BooleanFunction(N))
  L=[]
  while not wpb :
    E=Ekn(t,n)
    o=walsh_tr_0(f,E)//2
    #print(t)
    S=suppk(f,t)
    CS = [x for x in E if x not in S]
    L+=[(t,o,len(CS), binomial(len(CS),o))]
    for i in range(o):
      N[CS[i]]=1
    y=BooleanFunction(N)
    wpb, t =is_WPB(y) 
  nly=y.nonlinearity()
  print(nly,wpb)
  #print(L)
  return nly,y



def upmagic_rand(n):
  """Construction 1 from [GM22c] applied with magic(n) as input (bit-flipping randomised)"""
  f=magic(n)
  FT=f.truth_table('int')
  N=list(FT)
  N[-1]=1
  wpb, t =is_WPB(BooleanFunction(N))
  flipped={ n: [0] }
  while not wpb :
    E=Ekn(t,n)
    o=walsh_tr_0(f,E)//2
    #print(t)
    S=suppk(f,t)
    CS = [x for x in E if x not in S]
    O=sample(CS, o)
    flipped[t]=O
    for i in O:
      N[i]=1
    y=BooleanFunction(N)
    wpb, t =is_WPB(y) 
  nly=y.nonlinearity()
  print(nly,wpb)
  #if nly>=116: print(flipped)
  return nly, flipped, y
  
  
  
def collection_GM22c(n,d,outfile,B=0):
    """Creates a file including those functions obtained by applying upmagic_rand(n) d time whose
    NL is B. For n=8 and 16 B is set automatically to 116 and 32596, respectively.

    Args:
        n (int): number of variables
        d (int): number of calls of upmagic_rand(n)
        outfile (string): name of output file
    """    
    if n==8: B=116
    elif n==16: B=32596
    f = open(outfile, 'w')
    print('Output file:', f)
    print(outfile)
    f.write("---Constr 1---\n")
    f.write("n: "+str(n)+'\n')     
    stat={ 0: 0 }
    for i in range(d):
        nly, flipped, y =upmagic_rand(n)
        if nly in stat: stat[nly]+=1
        else: stat[nly]=1
        if nly>=B: 
          print(flipped)
          f.write(str(nly)+':'+str(flipped)+'\n')
    f.close()
    
    


  
