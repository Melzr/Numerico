# -*- coding: utf-8 -*-
import numpy as np


def newton_raphson(x0, f, df, max_iteraciones, tolerancia):

  """
  Calcula la raiz de la funcion f.
  Recibe una semilla x0, una funcion f y su derivada df continuas, una cantidad maxima de iteraciones y una tolerancia
  que sera la cota de la diferencia entre las dos ultimas iteraciones en caso de converger.
  Devolvera el historial de iteraciones con xn y el error como la diferencia con xn-1 en cada iteracion.
  Si la derivada de la funcion se anula, devolvera el resultado de la ultima iteracion.
  """

  if (df(x0) == 0):
    print("Error: division por cero")
    return None

  historial = []
  xn = (x0 - f(x0)/df(x0))
  xa = x0
  error = abs(xn - xa)

  #print("  n  |        Xn-1       |         Xn        |    Error")
  #print("    0|          -        | %.15f |     -" % (x0) )
  #print("    1| %.15f | %.15f | %.15f" % (x0, xn, error) )
  historial.append({'xn': xn, 'error': error})

  i = 2
  while ( (error > tolerancia) and (i < max_iteraciones) ):
    if (df(xn) == 0):
      print("Error: division por cero")
      break
    xa = xn
    xn = (xa - f(xa)/df(xa))
    error = abs(xn - xa)
    #print(" %4i| %.15f | %.15f | %.15f" % (i, xa, xn, error) )
    historial.append({'xn': xn, 'error': error})
    i+=1

  return historial


def biseccion_aux(px, a, b, f, tolerancia, iteraciones_restantes, historial):

  py = (a+b)/2

  if( abs(px - py) < tolerancia or iteraciones_restantes <= 0 ):
    return py
  
  if ( (f(py) > 0 and f(a) < 0) or (f(py) < 0 and f(a) > 0) ):
      b = py
  elif ( (f(py) > 0 and f(b) < 0) or (f(py) < 0 and f(b) > 0) ):
      a = py

  historial.append({'xn': py, 'error': abs(px - py)})
  return biseccion_aux(py, a, b, f, tolerancia, iteraciones_restantes-1, historial)


def biseccion(a, b, f, tolerancia, max_iteraciones):

  """
  Calcula la raiz de la funcion f.
  Recibe un intervalo [a,b], una funcion f, una tolerancia que sera la cota de la diferencia entre las dos ultimas
  iteraciones en caso de converger y una cantidad maxima de iteraciones.
  Devolvera el historial de iteraciones con xn y el error como la diferencia con xn-1 en cada iteracion.
  """

  if (f(a)*f(b) > 0):
    return None

  historial = []

  p0 = (a + b)/2
  historial.append({'xn': p0, 'error': abs(p0 - a)})
  if ( (f(p0) > 0 and f(a) < 0) or (f(p0) < 0 and f(a) > 0)):
    b = p0
  elif ( (f(p0) > 0 and f(b) < 0) or (f(p0) < 0 and f(b) > 0) ):
    a = p0
  
  biseccion_aux(p0, a, b, f, tolerancia, max_iteraciones, historial)
  return historial


def newton_raphson_modificado(x0, f, df, ddf, max_iteraciones, tolerancia):

  """
  Modificacion del metodo de Newton Raphson para funciones con raices multiples.
  A diferencia de Newton debe recibir ademas la derivada segunda de la funcion ddf.
  Devolvera el historial de iteraciones con xn y el error como la diferencia con xn-1 en cada iteracion.
  """

  if (df(x0) == 0):
    print("Error: division por cero")

  historial = []

  xn = (x0 - (f(x0) * df(x0))/(df(x0)**2 - f(x0)*ddf(x0)))
  xa = x0
  error = abs(xn - xa)
  historial.append({'xn': xn, 'error': error})
  #print("  n  |        Xn-1       |         Xn        |    Error")
  #print("    0|          -        | %.15f |     -" % (x0) )
  #print("    1| %.15f | %.15f | %.15f" % (x0, xn, error) )

  i = 2
  while ( (error > tolerancia) and (i < max_iteraciones) ):
    if (df(xn) == 0):
      print("Error: division por cero")
      break
    xa = xn
    xn = (xa - (f(xa) * df(xa))/(df(xa)**2 - f(xa)*ddf(xa)))
    error = abs(xn - xa)
    historial.append({'xn': xn, 'error': error})
    #print(" %4i| %.15f | %.15f | %.15f" % (i, xa, xn, error) )
    i+=1

  return historial


def secante(f, x0, x1, tolerancia, max_iteraciones):

  """
  Dada una funcion f continua, dos semillas x0 y x1, calcula su raiz de modo tal que la diferencia entre
  las dos ultimas iteraciones sera menor a la tolerancia.
  Devolvera el historial de iteraciones con xn y el error como la diferencia con xn-1 en cada iteracion.
  """

    historial = []
    error = abs(x0-x1)

    i = 0
    pn = x1
    historial.append({'xn': pn, 'error': error})

    while error > tolerancia and i < max_iteraciones :
        if (f(x1) == f(x0)):
            print("Error: division por cero")
            break
        pn = x1 - ( (f(x1)*(x1-x0))/(f(x1)-f(x0)) )
        x0=x1
        x1=pn
        error = abs(x0-x1)
        historial.append({'xn': pn, 'error': error})
    
    return historial
