#!/usr/bin/python3
## SageMath version 9.5

def randint_par(n):
  return np.random.RandomState().randint(0,n)


def partition(n):
    P={}
    for i in range(n+1): P[i]=[]
    for d in range(2^n): P[bin(d).count("1")].append(d)
    return P
    
def ver_p(n):
    P=partition(n)
    for k in range(n): 
      if len(P[k])!=binomial(n,k): return False
    print(True)

def it_bal_supp(b):
  return ithp(b,b//2)

def numWPB(m):
  n=2^m
  p=1
  for k in range(1,n): 
    b=binomial(n,k)
    p=p*binomial(b,b//2)
  return p

def WPB_iter(m):
  n=2^m
  Lit=[it_bal_supp(binomial(n,k)) for k in range(1,n)]
  return itertools.product(*Lit)

def wpb_from_it(m,S,P):
    n=2^m
    LF=2^n*[0]
    LF[-1]=1
    for k in range(1,n):
      suppk=S[k-1]
      indexk=P[k]
      for j in suppk: LF[indexk[j]]=1
    return  BooleanFunction(LF)
    
    
def stat_AI(m):
  P=partition(2^m)
  W=WPB_iter(m)
  n=2^m
  AI=n*[0]
  NL=2^(n-1)*[0]
  D=n*[0]
  for s in W: 
    Fs=wpb_from_it(m,s,P)
    assert is_WPB(Fs)
    AI[my_algebraic_immunity(Fs)]+=1
    NL[Fs.nonlinearity()]+=1
    D[Fs.algebraic_normal_form().degree()]+=1
  return AI,NL,D

def rand_WPB(m,P):
  n=2^m
  LF=2^n*[0]
  LF[-1]=1
  for k in range(1,n):
    I=set()
    indexk=P[k]
    b=binomial(n,k)
    while len(I)<b//2:
      j=randint_par(b)
      I.add(j)
      LF[indexk[j]]=1  
  return BooleanFunction(LF)
  
def stat_AI_rand(m,ns=1,verbose=True):
  P=partition(2^m)
  n=2^m
  AI=(n//2+1)*[0]
  NL=2^(n-1)*[0]
  D=n*[0]
  #Fs0=0
  for s in range(ns): 
    Fs=rand_WPB(m,P)
    #assert Fs0!=Fs
    #Fs0=Fs
    print(Fs)
    nl,d=Fs.nonlinearity(),Fs.algebraic_normal_form().degree()
    if verbose: print("NL : ", nl, "DEG : ",d, end=", ")
    a=my_algebraic_immunity(Fs)
    if verbose: print("AI : ", a)
    AI[a]+=1
    NL[nl]+=1
    D[d]+=1
  return AI,NL,D
  
  
#loop for the parallelisation of the random distribution estimation, see definition of test_rand_para 
def loop_rand_WPB(m,P,I): 
  Fs=rand_WPB(m,P)
  assert is_WPB(Fs)
  #print(Fs.truth_table('int'))
  #a=my_algebraic_immunity(Fs), Fs.nonlinearity(),Fs.algebraic_normal_form().degree()
  a=0,0,Fs.algebraic_normal_form().degree()
  print(a, end=',')
  return a
   


#use for print 'n_v_'+now_str()+'.txt'
def distribution_rand_para_all(m,n_sample=1000,coeff_par=cpu,verbose=False,outfile=None):
    print(cpu)
    P=partition(2^m)
    n=2^m
    AI=(n//2+1)*[0]
    NL=2^(n-1)*[0]
    D=n*[0]
    t0=time.time()
    treasure=[]
    i=0
    n_iter=int(n_sample/coeff_par)
    if outfile!=None:
      f = open(outfile, 'w')
      print('Output file:', f)
      print(outfile)
      f.write("---Random sample distribution---\n")
      f.write("n: "+str(n)+",m :"+str(m)+",n_sample: "+str(n_iter*coeff_par)+'\n')     
    while i<n_iter:
        print(i, end=',')
        pool = Pool(cpu)
        fo=partial(loop_rand_WPB,m,P)
        LV=pool.map(fo,range(coeff_par))
        #print(LV)
        pool.close()
        pool.join()
        for c in LV: 
          AI[c[0]]+=1
          NL[c[1]]+=1
          D[c[2]]+=1
        del LV
        i+=1
    t1= time.time()   
    print(i)
    print("\nrunning time: %.2f" %(t1-t0))
    
    if outfile!=None:
        f.write("\n\nAI: "+str(AI))
        f.write("\n\nNL: "+str(NL))
        f.write("\n\nDEG: "+str(D))
        f.write("\nrunning time: %.2f sec" %(t1-t0))
        f.close()
    return vector(ZZ,AI),vector(ZZ,NL),vector(ZZ,D) 

"""
m=3
AI    (0, 0, 0, 2, 9982, 0, 0, 0)
NL (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 9, 0, 21, 0, 58, 0, 167, 0, 454, 0, 1004, 0, 2039, 0, 2967, 0, 2510, 0, 723, 0, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

"""
