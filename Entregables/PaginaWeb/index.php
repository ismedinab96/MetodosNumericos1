<!DOCTYPE html>
<html lang="es">
<head>
  <meta charset="UTF-8">
  <title>Simulador Fotovoltaico</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
  <link href="https://cdn.jsdelivr.net/npm/bootstrap-icons@1.11.1/font/bootstrap-icons.css" rel="stylesheet">

  <style>
    body {
      background: linear-gradient(135deg, #f8f9fa, #e3f2fd);
    }
    .hero {
      background: url('https://source.unsplash.com/1600x400/?solar,energy') center/cover no-repeat;
      color: white;
      padding: 4rem 2rem;
      text-shadow: 2px 2px 4px rgba(0,0,0,0.6);
    }
    .card {
      border-radius: 12px;
      box-shadow: 0 4px 12px rgba(0,0,0,0.1);
    }
  </style>
</head>

<body>

<div class="hero text-center">
  <h1 class="display-4 fw-bold">⚡ Metodos Numericos I</h1>
  <h1 class="display-4 fw-bold">⚡ Microrred de Sistemas Fotovoltaicos</h1>
  <p class="lead">Calculadora para resolver sistemas de ecuaciones lineales</p>
</div>

<div class="container mt-5">

  <!-- FORMULARIO -->
  <div class="card p-4">
    <h3 class="text-center mb-4">Ingrese los coeficientes del sistema:</h3>

    <form method="POST">

      <div class="row">

        <!-- MATRIZ A -->
        <div class="col-md-8">
          <h5 class="text-center">Matriz A</h5>
          <table class="table table-bordered text-center">
            <tbody>
              <?php
              for($i=0;$i<5;$i++){
                echo "<tr>";
                for($j=0;$j<5;$j++){
                  echo "<td>
                    <input class='form-control text-center'
                           type='number' step='any'
                           name='A$i$j' required>
                  </td>";
                }
                echo "</tr>";
              }
              ?>
            </tbody>
          </table>
        </div>

        <!-- VECTOR b -->
        <div class="col-md-4">
          <h5 class="text-center">Vector b</h5>
          <table class="table table-bordered text-center">
            <tbody>
              <?php
              for($i=0;$i<5;$i++){
                echo "<tr>
                  <td>
                    <input class='form-control text-center'
                           type='number' step='any'
                           name='b$i' required>
                  </td>
                </tr>";
              }
              ?>
            </tbody>
          </table>
        </div>

      </div>

      <div class="text-center mt-4">
        <button type="submit" class="btn btn-success btn-lg">
          <i class="bi bi-lightning-charge-fill"></i> Resolver Sistema
        </button>
      </div>

    </form>
  </div>

  <!-- PHP -->
  <?php
  error_reporting(E_ALL);
  ini_set('display_errors', 1);

  if($_SERVER["REQUEST_METHOD"] == "POST"){

    $A = [];
    $b = [];

    // MATRIZ A
    for($i=0;$i<5;$i++){
      for($j=0;$j<5;$j++){
        $key = "A".$i.$j;
        if(!isset($_POST[$key]) || $_POST[$key] === ""){
          echo "<div class='alert alert-danger mt-4'>Falta dato en A[$i][$j]</div>";
          return;
        }
        $A[$i][$j] = floatval($_POST[$key]);
      }
    }

    // VECTOR b
    for($i=0;$i<5;$i++){
      $key = "b".$i;
      if(!isset($_POST[$key]) || $_POST[$key] === ""){
        echo "<div class='alert alert-danger mt-4'>Falta dato en b[$i]</div>";
        return;
      }
      $b[$i] = floatval($_POST[$key]);
    }

    // Validación diagonal
    for($i=0;$i<5;$i++){
      if($A[$i][$i] == 0){
        echo "<div class='alert alert-danger mt-4'>División por cero en fila ".($i+1)."</div>";
        return;
      }
    }

    // -------- Gauss-Seidel --------
function gs($A,$b,$tol=1e-6,$maxIter=100){
  $x=array_fill(0,5,0);
  $iter=0;

  do{
    $x_old=$x;

    for($i=0;$i<5;$i++){
      $sum=0;
      for($j=0;$j<5;$j++){
        if($j!=$i) $sum+=$A[$i][$j]*$x[$j];
      }
      $x[$i]=($b[$i]-$sum)/$A[$i][$i];
    }

    $error=0;
    for($i=0;$i<5;$i++){
      $error=max($error, abs($x[$i]-$x_old[$i]));
    }

    $iter++;

  }while($error>$tol && $iter<$maxIter);

  return ["x"=>$x,"iter"=>$iter,"error"=>$error];
}

    // -------- SOR --------
function sor($A,$b,$w=1.1,$tol=1e-6,$maxIter=100){
  $x=array_fill(0,5,0);
  $iter=0;

  do{
    $x_old=$x;

    for($i=0;$i<5;$i++){
      $sum=0;
      for($j=0;$j<5;$j++){
        if($j!=$i) $sum+=$A[$i][$j]*$x[$j];
      }
      $x[$i]=(1-$w)*$x[$i]+$w*(($b[$i]-$sum)/$A[$i][$i]);
    }

    $error=0;
    for($i=0;$i<5;$i++){
      $error=max($error, abs($x[$i]-$x_old[$i]));
    }

    $iter++;

  }while($error>$tol && $iter<$maxIter);

  return ["x"=>$x,"iter"=>$iter,"error"=>$error];
}

    // -------- Gradiente Conjugado --------
function cg($A,$b,$tol=1e-6,$maxIter=100){
  $x=array_fill(0,5,0);

  for($i=0;$i<5;$i++){
    $Ax=0;
    for($j=0;$j<5;$j++) $Ax+=$A[$i][$j]*$x[$j];
    $r[$i]=$b[$i]-$Ax;
    $p[$i]=$r[$i];
  }

  $iter=0;

  do{
    for($i=0;$i<5;$i++){
      $Ap[$i]=0;
      for($j=0;$j<5;$j++) $Ap[$i]+=$A[$i][$j]*$p[$j];
    }

    $rTr=0; $pAp=0;
    for($i=0;$i<5;$i++){
      $rTr+=$r[$i]*$r[$i];
      $pAp+=$p[$i]*$Ap[$i];
    }

    if($pAp==0) break;

    $alpha=$rTr/$pAp;

    for($i=0;$i<5;$i++){
      $x[$i]+=$alpha*$p[$i];
      $r_new[$i]=$r[$i]-$alpha*$Ap[$i];
    }

    $error=0;
    for($i=0;$i<5;$i++){
      $error=max($error, abs($r_new[$i]));
    }

    if($error<$tol) break;

    $rTr_new=0;
    for($i=0;$i<5;$i++) $rTr_new+=$r_new[$i]*$r_new[$i];

    $beta=$rTr_new/$rTr;

    for($i=0;$i<5;$i++){
      $p[$i]=$r_new[$i]+$beta*$p[$i];
      $r[$i]=$r_new[$i];
    }

    $iter++;

  }while($iter<$maxIter);

  return ["x"=>$x,"iter"=>$iter,"error"=>$error];
}

$gs = gs($A,$b);
$sor = sor($A,$b);
$cg = cg($A,$b);

function mostrar($res){
  $out="";
  foreach($res["x"] as $i=>$v){
    $out.="x".($i+1)." = ".round($v,4)."<br>";
  }
  $out.="<hr>Iteraciones: ".$res["iter"];
  $out.="<br>Error: ".round($res["error"],8);
  return $out;
}

    echo '
    <div class="card mt-5 shadow-lg">
      <div class="card-body">
        <h4 class="mb-3 text-primary">Resultados</h4>

        <ul class="nav nav-tabs">
          <li class="nav-item">
            <button class="nav-link active" data-bs-toggle="tab" data-bs-target="#gs">Gauss-Seidel</button>
          </li>
          <li class="nav-item">
            <button class="nav-link" data-bs-toggle="tab" data-bs-target="#sor">SOR</button>
          </li>
          <li class="nav-item">
            <button class="nav-link" data-bs-toggle="tab" data-bs-target="#cg">Gradiente Conjugado</button>
          </li>
        </ul>

        <div class="tab-content mt-3">
            <div class="tab-pane fade show active" id="gs">'.mostrar($gs).'</div>
            <div class="tab-pane fade" id="sor">'.mostrar($sor).'</div>
            <div class="tab-pane fade" id="cg">'.mostrar($cg).'</div>
        </div>

      </div>
    </div>';
  }
  ?>

</div>

<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>

<footer class="bg-dark text-center text-white mt-5">
  <div class="p-3" style="background-color: rgba(0, 0, 0, 0.2);">
    © 2026 Elaborado por: Medina Balboa Iver Samuel<br>
    <b>Métodos Numéricos I - Lic. Brigida Carvajal Blanco</b>
  </div>
</footer>

</body>
</html>