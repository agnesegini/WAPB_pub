import random

def elem_2(n):
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  f=R(0)
  for i in range(n):
    for j in range(i):
     f+=R('x%d*x%d' %(i,j))
  return f
      
def magic(n):
  R=BooleanPolynomialRing(names=['x'+str(i) for i in range(n)])
  f=sum([R('x%d' %i) for i in range(0,n//2)])
  for i in range(n):
    for j in range(i):
     f+=R('x%d*x%d' %(i,j))
  return BooleanFunction(f)

def comp(f,n):
  return sum([abs(walsh_tr_0(f,Ekn(i,n))) for i in range(n+1)])

def countm(n):
   f=magic(n)
   L=[walsh_tr_0(f,Ekn(i,n)) for i in range(n+1)]
   print(L)
   t=1
   for k in range(1,n):
      t=t*binomial(binomial(n,k)/2+L[k]/2,L[k])
   print(t,log(t,2).n())
   return t

def upmagic(n):
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
  return nly#,y



def upmagic_rand(n):
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
    O=random.sample(CS, o)
    flipped[t]=O
    for i in O:
      N[i]=1
    y=BooleanFunction(N)
    wpb, t =is_WPB(y) 
  nly=y.nonlinearity()
  print(nly,wpb)
  #if nly>=116: print(flipped)
  return nly, flipped#, y
  
  
  
def collection(n,d,outfile):
    if n==8: B=116
    elif n==16: B=32596
    f = open(outfile, 'w')
    print('Output file:', f)
    print(outfile)
    f.write("---Constr 1---\n")
    f.write("n: "+str(n)+'\n')     
    stat={ 0: 0 }
    for i in range(d):
        nly, flipped=upmagic_rand(n)
        if nly in stat: stat[nly]+=1
        else: stat[nly]=1
        if nly>=B: 
          print(flipped)
          f.write(str(nly)+':'+str(flipped)+'\n')
    f.close()
    
    
#collection(16,50,'mer-new.txt')

  
