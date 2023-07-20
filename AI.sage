def my_algebraic_immunity(f, annihilator = False, verbose = False):
        g = ~f
        n=f.nvariables()
        for i in range(ceil(n/2)):  #see carlet 9.1
            if verbose : print(i, end=': ')
            for fun in [f, g]:
                if verbose : print(fun, end=' - ')
                A = fun.annihilator(i)
                if A is not None:
                    if annihilator:
                        return i, A
                    else:
                        return i
                if verbose :  print("")
        if verbose : print(ceil(n/2))
        return(ceil(n/2)) 
        raise ValueError("you just found a bug!")
