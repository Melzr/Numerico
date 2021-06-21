# -*- coding: utf-8 -*-
import numpy as np


def orden_de_convergencia(pn, pn_mas_uno, pn_menos_uno, pn_menos_dos):
    
  """
  Calcula el orden de convergencia de un metodo en una iteracion pn, las dos anteriores y la siguiente
  """
  
  error_pn = pn - pn_menos_uno
  error_pn_mas_uno = pn_mas_uno - pn
  error_pn_menos_uno = pn_menos_uno - pn_menos_dos

  if abs(error_pn) < 1e-14 or abs(error_pn_mas_uno) < 1e-14 or abs(error_pn_menos_uno) < 1e-14 or abs(np.log10(abs(error_pn/error_pn_menos_uno))) < 1e-14:
    return None

  return np.log10(abs(error_pn_mas_uno/error_pn))/np.log10(abs(error_pn/error_pn_menos_uno))


def constante_asistotica(pn, pn_mas_uno, p, orden_de_convergencia):

  """
  Calcula la constante de error asintotica de un metodo en una iteracion n siendo p la raiz de la funcion o una
  aproximacion de la misma.
  El orden de convergencia debe ser el de la misma iteracion n.
  """

  dif_pn = pn - pn_menos_uno
  if_pn_mas_uno = pn_mas_uno - pn 
  if  orden_de_convergencia is None or (abs(dif_pn)**orden_de_convergencia) < 1e-14 or abs(dif_pn_mas_uno)/(abs(dif_pn)**orden_de_convergencia) > 1 or abs(dif_pn_mas_uno)/(abs(dif_pn)**orden_de_convergencia) < 0.1:
    return None
        
  return abs(dif_pn_mas_uno)/(abs(dif_pn)**orden_de_convergencia)


def calcular_constante_asintotica(historial):

  """
  Dado el historial de un metodo de calculo de raiz con el valor de la raiz y su error en cada iteracion
  devuelve el historial de la constante asintotica en cada iteracion.
  """

  ns = [None, None, None]
  while len(historial) > 4:
    convergencia = orden_de_convergencia((historial[2])['xn'], (historial[3])['xn'], (historial[1])['xn'], (historial[0])['xn'])
    ns.append(constante_asistotica((historial[2])['xn'],(historial[3])['xn'],(historial[1])['xn'], (historial[-1])['xn'], convergencia))
    historial.pop(0)
  return ns


def calcular_convergencia(historial):

  """
  Dado el historial de un metodo de calculo de raiz con el valor de la raiz y su error en cada iteracion
  devuelve el historial del orden de convergencia en cada iteracion.
  """

  ns = [None, None, None]
  while len(historial) > 4:
    convergencia = orden_de_convergencia((historial[2])['xn'], (historial[3])['xn'], (historial[1])['xn'], (historial[0])['xn'])
    ns.append(convergencia)
    historial.pop(0)
  return ns