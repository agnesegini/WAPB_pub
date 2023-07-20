#!/usr/bin/python3
## SageMath version 9.5
#load("main.sage")
import time
import numpy as np

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

def it_bal_supp_uptosymm(b):
  return ithp(b-1,b//2-1)
  


def numWPB(m):
  n=2^m
  p=1
  for k in range(1,n): 
    b=binomial(n,k)
    p=p*binomial(b,b//2)
  return p

def numWPB_uptosymm(m):
  n=2^m
  p=1
  for k in range(1,n): 
    b=binomial(n,k)
    p=p*binomial(b,b//2)
  return p/(2^(n-1))

def WPB_iter_uptosymm(m):
  n=2^m
  Lit=[it_bal_supp_uptosymm(binomial(n,k)) for k in range(1,n)]
  return itertools.product(*Lit)

#canonical rapresentation with last element of the slice having always true evaluation
def wpb_from_it_uptosymm(m,S,P): 
    n=2^m
    LF=2^n*[0]
    LF[-1]=1
    for k in range(1,n):
      suppk=S[k-1]
      indexk=P[k]
      for j in suppk: LF[indexk[j]]=1
      LF[indexk[-1]]=1
    return  BooleanFunction(LF)
    
def rand_WPB_uptosymm(m,P=0):
  n=2^m
  if P==0: P=partition(n)
  LF=2^n*[0]
  LF[-1]=1
  for k in range(1,n):
    I=set()
    indexk=P[k]
    b=binomial(n,k)
    while len(I)<b//2-1:
      j=randint_par(b-1)
      I.add(j)
      LF[indexk[j]]=1  
    LF[indexk[-1]]=1
  return BooleanFunction(LF)
  
def sub_to_symm(x,I,i,f):
   for j in range(i): 
      if x[j]: f+=I[j]
   return f
   
def s_class_0(n,P=0):
   if P==0: P=partition(n)
   I=[Ikn(k,n) for k in range(1,n)]
   f=BooleanFunction(list(zero_vector(2^n))) 
   i=len(I)
   sub_it=itertools.product(range(2),repeat=i)
   return map(lambda x : sub_to_symm(x,I,i,f), sub_it)
      
def s_class(f,P=0):
   n=f.nvariables()
   if P==0: P=partition(n)
   return [f+s for s in s_class_0(n,P)]

def stat_S_class(f,S0=0,full=True):
    n=f.nvariables()
    if S0==0: S0 = s_class_0(n)
    L=[]
    NL,AI,DEG=set(),set(),set()
    for s in S0: 
      g=f+s
      #assert is_WPB(g)
      nl,ai,d = g.nonlinearity(),g.algebraic_immunity(),g.algebraic_normal_form().degree()
      NL.add(nl)
      AI.add(ai)
      DEG.add(d)
      #print(nl,ai,d)
      if full: L+=[(nl,ai,d)]
      del g
    print(NL,AI,DEG)
    return f, NL,AI,DEG,L
    
def stat_S_class_all(m=2):
  P=partition(2^m)
  n=2^m
  S0 = list(s_class_0(n))
  L=[]
  for U in WPB_iter_uptosymm(m):
   f=wpb_from_it_uptosymm(m,U,P)
   Sf=stat_S_class(f,S0,False)
   rf=[NLk(k,f) for k in range(1,n)]
   print(Sf, rf)
   L+=[[Sf,rf]]
  return L
  

def taxonomy4(table=True):
   L=stat_S_class_all(2)
   print("!")
   T=dict()
   pr=[]
   for Sf in L: 
      a=str([Sf[1][1],Sf[0][1],Sf[0][3]])
      if a in pr: T[a]+=1
      else: 
         T[a]=1
         pr.append(a)
   if table:
      for a in pr: print( a, ' & ',  T[a], " \\\\")  
         
def stat_S_class_rand(m,ns=1):
  P=partition(2^m)
  n=2^m
  S0 = list(s_class_0(n))
  L=[]
  for s in range(ns): 
    f=rand_WPB_uptosymm(m,P)
    L+=[stat_S_class(f,P,S0,False)]
  return L

  
#loop for the parallelisation of
def loop_rand_S_class(Fun,s): 
  #print(Fs.truth_table('int'))
  Fs=Fun+s
  a=my_algebraic_immunity(Fs), Fs.nonlinearity(),Fs.algebraic_normal_form().degree()
  #print(a, end=',')
  return a
   


#use for print 'n_v_'+now_str()+'.txt'
def stat_S_class_rand_para(m,verbose=False,outfile=None):
    print(cpu)
    P=partition(2^m)
    n=2^m
    Fun=rand_WPB_uptosymm(m,P)
    NL,AI,DEG=set(),set(),set()
    distAI=(n//2+2)*[0]
    distNL=2^(n-1)*[0]
    distD=n*[0]
    t0=time.time()
    treasure=[]
    if outfile!=None:
      f = open(outfile, 'w')
      print('Output file:', f)
      print(outfile)
      f.write("---S-Class---\n")
      f.write("n: "+str(n)+",m :"+str(m)+'\n')     
    pool = Pool(cpu)
    fo=partial(loop_rand_S_class,Fun)
    S0=s_class_0(n)
    LV=pool.map(fo,S0)
   #print(LV)
    pool.close()
    pool.join()
    for c in LV:
          NL.add(c[1])
          AI.add(c[0])
          DEG.add(c[2])
          distAI[c[0]]+=1
          distNL[c[1]]+=1
          distD[c[2]]+=1
    t1= time.time()   
    print("\nrunning time: %.2f" %(t1-t0))
    if outfile!=None:
        f.write("\n\ndistAI: "+str(distAI))
        f.write("\n\nAI: "+str(AI))
        f.write("\n\ndistNL: "+str(distNL))
        f.write("\n\nNL: "+str(NL))
        f.write("\n\ndistDEG: "+str(distD))
        f.write("\n\ndistDEG: "+str(DEG))
        f.write("\nclass representative: " + str(Fun.truth_table(format='int')))
        f.close()
    return vector(ZZ,distAI),vector(ZZ,distNL),vector(ZZ,distD),AI,NL,DEG 


def stat_S_class_para(Fun,verbose=False,outfile=None):
    print(cpu)
    n=Fun.nvariables()
    NL,AI,DEG=set(),set(),set()
    distAI=(n//2+2)*[0]
    distNL=2^(n-1)*[0]
    distD=n*[0]
    t0=time.time()
    treasure=[]
    if outfile!=None:
      f = open(outfile, 'w')
      print('Output file:', f)
      print(outfile)
      f.write("---S-Class---\n")
      f.write("n: "+str(n)+'\n')     
    pool = Pool(cpu)
    fo=partial(loop_rand_S_class,Fun)
    S0=s_class_0(n)
    LV=pool.map(fo,S0)
   #print(LV)
    pool.close()
    pool.join()
    for c in LV:
          NL.add(c[1])
          AI.add(c[0])
          DEG.add(c[2])
          distAI[c[0]]+=1
          distNL[c[1]]+=1
          distD[c[2]]+=1
    t1= time.time()   
    print("\nrunning time: %.2f" %(t1-t0))
    if outfile!=None:
        f.write("\n\ndistAI: "+str(distAI))
        f.write("\n\nAI: "+str(AI))
        f.write("\n\ndistNL: "+str(distNL))
        f.write("\n\nNL: "+str(NL))
        f.write("\n\ndistDEG: "+str(distD))
        f.write("\n\nDEG: "+str(DEG))
        
        f.write("\nrunning time: %.2f sec" %(t1-t0))
        
        f.write("\nclass representative: " + str(Fun.truth_table(format='int')))
        f.close()
    return vector(ZZ,distAI),vector(ZZ,distNL),vector(ZZ,distD),AI,NL,DEG 
