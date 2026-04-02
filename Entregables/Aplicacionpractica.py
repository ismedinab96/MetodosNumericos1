import numpy as np
#Ecuaciones lineales en su forma matricial de la forma Ax=b
A = np.array([[10,-2,-1,0,0],[-2,10,-2,-1,0],[-1,-2,10,-2,-1],[0,-1,-2,10,-2],[0,0,-1,-2,10]], dtype = float)
b = np.array([7,8,6,5,4], dtype=float)
#Valor inicial de ceros para el metodo
x0 = np.zeros(5)

# Aplicacion del Metodo de Gauss Seidel

def met_gauss_Seidel(A,b,x, it=15):
    n = len(b)
    x=x.copy()
    
    print(f'Metodo de Gauss Seidel con {it} iteraciones')
    for k in range(it):
        for i in range(n):
            s1=sum(A[i,j]*x[j] for j in range(i))
            s2=sum(A[i,j]*x[j] for j in range(i+1,n))
            
            x[i]=(b[i]-s1-s2)/A[i,i]
        print(f'Iteraciones: {k+1} -> {x}')

met_gauss_Seidel(A,b,x0)

# Aplicacion del Metodo de SOR

def met_SOR(A,b,x,w=1.1,it=15):
    n = len(b)
    x = x.copy()
    
    print(f'Metodo de SOR con {it} iteraciones y omega = {w}')
    
    for k in range(it):
        for i in range(n):
            s1=sum(A[i,j]*x[j] for j in range(i))
            s2=sum(A[i,j]*x[j] for j in range (i+1,n))
            
            xi=(b[i]-s1-s2)/A[i,i]
            
            x[i]=(1-w)*x[i] + w*xi
        print(f'Iteraciones: {k+1} -> {x}')
met_SOR(A,b,x0)

# Aplicacion del metodo de gradiente conjugado

def met_gradiente_conjugado(A,b,x,it=15):
    x = x.copy()
    
    r= b - A@x
    p = r.copy()
    
    print(f'Metodo de Gradiente Conjugado con {it} Iteraciones')
    for k in range(it):
        Ap = A @ p
        alpha = (r @ r) / (p @ Ap)
        x = x + alpha*p
        r2 = r - alpha*Ap
        print(f'Iteraciones : {k+1} -> {x}')
        
        if np.linalg.norm(r2) < 1e-10:
            break
        beta = (r2 @ r2) / (r @ r)
        p = r2 + beta*p
        r = r2
    return x

met_gradiente_conjugado(A,b,x0)