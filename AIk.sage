

def build_G(n,d):
        P = partition(n)
        C = codes.BinaryReedMullerCode(d, n)
        G = C.generator_matrix()
        GK= [G[:, P[i]] for i in range(0,n+1) ]
        return GK

def weightwise_annihilator(f,k,d,G=0,P=0):
    n = f.nvariables()
    if not P: P = partition(n)
    if not G: G=build_G(n,d)
    Gk,Pk=G[k],P[k]
    ke=Gk.left_kernel().dimension()
    M0,M1=[],[]
    for i in range(len(P[k])):
          if f(Pk[i]): M1+=[Gk.T[i]]
          else: M0+=[Gk.T[i]]
    M0,M1=matrix(M0).T,matrix(M1).T
    A0,A1=M0.kernel(),M1.kernel()
    a0,a1=A0.dimension()-ke,A1.dimension()-ke
    return a0,a1, A0,A1
    
def AIk(k,f):
    n = f.nvariables()
    P = partition(n)
    a0,a1,d =0,0,1
    while a0==0 and a1==0 and d<n:
        a0,a1,A0,A1=weightwise_annihilator(f,k,d,G=0,P=P)
        print(a0,a1)
        d+=1
    return d-1

def AIk_all(f):
    n = f.nvariables()
    P = partition(n)
    A=zero_vector(n+1)
    d=1
    L=list(range(1,n))
    while L!=[]:
        L1=[]
        G=build_G(n,d) 
        for k in L: 
             a0,a1,A0,A1=weightwise_annihilator(f,k,d,G,P)
             if a0>0 or a1>0: 
                  A[k]=d
             else: L1.append(k)
        L=L1
        d+=1
    return A

             
     
"""
 def build_G_2(n,d):
        P = partition(n)
        C = codes.BinaryReedMullerCode(d, n)
        GK= []
        R=range(2**n)
        for i in range(0,n+1):
              S=set(R).difference(set(P[i]))
              Cp=codes.PuncturedCode(C, S)
              #print(Cp.dimension(),Cp.length())
              GK+=[Cp.generator_matrix()]
        return GK   
"""

