from flask import Flask, render_template, request

app = Flask(__name__)

# --------- MÉTODOS ---------

def gauss_seidel(A, b, it=50):
    x = [0]*len(b)
    for _ in range(it):
        for i in range(len(b)):
            s = sum(A[i][j]*x[j] for j in range(len(b)) if j != i)
            x[i] = (b[i] - s)/A[i][i]
    return x

def sor(A, b, w=1.1, it=50):
    x = [0]*len(b)
    for _ in range(it):
        for i in range(len(b)):
            s = sum(A[i][j]*x[j] for j in range(len(b)) if j != i)
            x[i] = (1-w)*x[i] + w*((b[i] - s)/A[i][i])
    return x

def gradiente_conjugado(A, b):
    n = len(b)
    x = [0]*n

    def mult(A, x):
        return [sum(A[i][j]*x[j] for j in range(n)) for i in range(n)]

    r = [b[i] - mult(A,x)[i] for i in range(n)]
    p = r[:]

    for _ in range(n):
        Ap = mult(A,p)
        rTr = sum(r[i]*r[i] for i in range(n))
        pAp = sum(p[i]*Ap[i] for i in range(n))
        alpha = rTr / pAp

        x = [x[i] + alpha*p[i] for i in range(n)]
        r_new = [r[i] - alpha*Ap[i] for i in range(n)]

        rTr_new = sum(r_new[i]*r_new[i] for i in range(n))
        beta = rTr_new / rTr

        p = [r_new[i] + beta*p[i] for i in range(n)]
        r = r_new

    return x

# --------- RUTA WEB ---------

@app.route("/", methods=["GET","POST"])
def index():

    if request.method == "POST":

        A = []
        b = []

        for i in range(5):
            fila = []
            for j in range(5):
                fila.append(float(request.form[f"a{i}{j}"]))
            A.append(fila)
            b.append(float(request.form[f"a{i}5"]))

        gs = gauss_seidel(A,b)
        sor_res = sor(A,b)
        cg = gradiente_conjugado(A,b)

        return render_template("index.html", gs=gs, sor=sor_res, cg=cg)

    return render_template("index.html")

# IMPORTANTE PARA RENDER
if __name__ == "__main__":
    app.run()s