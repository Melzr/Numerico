# -*- coding: utf-8 -*-
import numpy as np


def orden_de_convergencia(error_pn, error_pn_mas_uno, error_pn_menos_uno):
    
    """
    Calcula el orden de convergencia de un metodo en una iteracion pn, dados los errores que son la diferencia
    de la iteracion correspondiente con la anterior.
    """

  if abs(error_pn) < 1e-14 or abs(error_pn_mas_uno) < 1e-14 or abs(error_pn_menos_uno) < 1e-14 or abs(np.log10(abs(error_pn/error_pn_menos_uno))) < 1e-14:
    return None

  return np.log10(abs(error_pn_mas_uno/error_pn))/np.log10(abs(error_pn/error_pn_menos_uno))


def constante_asistotica(pn, pn_mas_uno, p, orden_de_convergencia):

    """
    Calcula la constante de error asintotica de un metodo en una iteracion n siendo p la raiz de la funcion o una
    aproximacion de la misma.
    El orden de convergencia debe ser el de la misma iteracion n.
    """
    if pn == p or orden_de_convergencia is None:
        return None
        
    return abs(pn_mas_uno - p)/(abs(pn - p)**orden_de_convergencia)


def calcular_constante_asintotica(historial):

    """
    Dado el historial de un metodo de calculo de raiz con el valor de la raiz y su error en cada iteracion
    devuelve el historial de la constante asintotica en cada iteracion.
    """

  ns = [None, None]
  while len(historial) > 3:
    convergencia = orden_de_convergencia((historial[1])['error'], (historial[2])['error'], (historial[0])['error'])
    ns.append(constante_asistotica((historial[1])['xn'],(historial[2])['xn'], (historial[-1])['xn'], convergencia))
    historial.pop(0)
  return ns


def calcular_convergencia(historial):

    """
    Dado el historial de un metodo de calculo de raiz con el valor de la raiz y su error en cada iteracion
    devuelve el historial del orden de convergencia en cada iteracion.
    """

  ns = [None, None]
  while len(historial) > 3:
    convergencia = orden_de_convergencia((historial[1])['error'], (historial[2])['error'], (historial[0])['error'])
    ns.append(convergencia)
    historial.pop(0)
  return ns
  