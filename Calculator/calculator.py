#!/usr/bin/env python3

import sys,os
from string import *
from tkinter import *
from PIL import Image
from math import *
from time import *
from random import *
#from Canvas import *
from scipy.integrate import quad
from scipy.optimize import leastsq
from scipy.optimize import bisect
import numpy as np
import matplotlib.pyplot as plt

from tkinter import messagebox

def onclick1():
  pole.insert(END, "1")
def onclick2():
  pole.insert(END, "2")
def onclick3():
  pole.insert(END, "3")

def onclick4():
  pole.insert(END, "4")
def onclick5():
  pole.insert(END, "5")
def onclick6():
  pole.insert(END, "6")

def onclick7():
  pole.insert(END, "7")
def onclick8():
  pole.insert(END, "8")
def onclick9():
  pole.insert(END, "9")

def onclick0():
  pole.insert(END, "0")
def onclickpoint():
  pole.insert(END, ".")


def onclickpl():
  pole.insert(END, "+")
def onclickmin():
  pole.insert(END, "-")
def onclicktimes():
  pole.insert(END, "*")
def onclickdiv():
  pole.insert(END, "/")
def onclickpow():
  pole.insert(END, "**")
def onclickmod():
  pole.insert(END, "%")

def onclickdel():
  s=pole.get()
  s=s[0:-1]
  pole.delete(0,END)
  pole.insert(0,s)

def onclickrow():
  try:
    wart=eval(pole.get())
    pole.delete(0,END)
    pole.insert(0, str(wart))
  except:
    pole.delete(0,END)
    pole.insert(0, "ERROR")

def onclickCalka():
  s=pole.get()
  s=s.split()

  try:
    q=float(s[0])
    messagebox.showinfo("Info","Podaj poprawną nazwę funkcji")
  except ValueError:
    pass

  if len(s)!=3:
    messagebox.showinfo("Info","Podaj funkcję oraz granice")
    return 0


  try:
    a=eval(s[1])
    b=eval(s[2])
  except:
    messagebox.showinfo("Błąd","Granice muszą być wyrażone liczbą")
    return 0
  if a>=b:
    messagebox.showinfo("Błąd","Błędne granice")
    return 0

  fun=lambda x: eval(s[0])
  try:
    q=quad(fun,a,b)
    pole.delete(0,END)
    pole.insert(0,str(q[0]))

  except:
    messagebox.showinfo("Błąd","Wystąpił błąd")
    return 0



def onclickZero():
  s=pole.get()
  s=s.split()

  try:
    q=float(s[0])
    messagebox.showinfo("Info","Podaj poprawną nazwę funkcji")
  except ValueError:
    pass

  if len(s)!=3:
    messagebox.showinfo("Info","Podaj funkcję oraz granice")
    return 0


  try:
    a=eval(s[1])
    b=eval(s[2])
  except:
    messagebox.showinfo("Błąd","Granice muszą być wyrażone liczbą")
    return 0
  if a>=b:
    messagebox.showinfo("Błąd","Błędne granice")
    return 0

  fun=lambda x: eval(s[0])


  try:
    q=bisect(fun,a,b)
    pole.delete(0,END)
    pole.insert(0,str(q))

  except:
    messagebox.showinfo("Błąd","Podaj inne granice")
    return 0



def onclickLoad():
  global foto
  plt.clf()
  
  s=pole.get()
  s=s.split()

  try:
    q=float(s[0])
    messagebox.showinfo("Info","Podaj poprawną nazwę pliku")
  except ValueError:
    pass

  if len(s)!=1:
    messagebox.showinfo("Info","Podaj nazwę pliku")
    return 0

  try:
    X,Y=np.loadtxt(s[0],unpack=True)
    """
    X=[i[0] for i in f]
    Y=[i[1] for i in f]
    X=np.array(X)
    Y=np.array(Y)
    """
  except:
    messagebox.showinfo("Błąd","Nie udało się wczytać danych")
    return 0

  func=lambda t,x: t[0]*x**2+t[1]*x+t[2]
  errfunc=lambda t,x,y: func(t,x)-y
  ti=(1.0,2.0,0.0)
  tfin, suc = leastsq(errfunc,ti[:],args=(X,Y))

  xx=np.linspace(X.min(),X.max(),50)
  yy=func(tfin,xx)

  plt.plot(X,Y,'b*',xx,yy,'r')
  plt.savefig('wykr.png', dpi=65)
  s='a = '+str(tfin[0])+' b = '+str(tfin[1])+' c = '+str(tfin[2])
  
  try:
   Image.open('wykr.png').save('wykr.gif')
   foto=PhotoImage(file="wykr.gif")
   cv.create_image(0,0,anchor=NW,image=foto)
  except:
    pass
  messagebox.showinfo("Wynik dopasowania",s)

def onclickRysuj():
  global foto
  plt.clf()
  s=pole.get()
  s=s.split()

  try:
    q=float(s[0])
    messagebox.showinfo("Info","Podaj poprawną nazwę funkcji")
  except ValueError:
    pass

  if len(s)!=3:
    messagebox.showinfo("Info","Podaj funkcję oraz granice")
    return 0


  try:
    a=eval(s[1])
    b=eval(s[2])
  except:
    messagebox.showinfo("Błąd","Granice muszą być wyrażone liczbą")
    return 0
  if a>=b:
    messagebox.showinfo("Błąd","Błędne granice")
    return 0
  try:
    X=np.arange(a,b,0.01)
    #fun=lambda x: eval(s[0])
    Y=[eval(s[0]) for x in X]
    Y=np.array(Y)
    #plt.plot(X,fun(X))
    plt.plot(X,Y)
  except ValueError:
    messagebox.showinfo("Błąd","Błędne granice")


  plt.savefig('wykr.png', dpi=65)
  
  try:
   Image.open('wykr.png').save('wykr.gif')
   foto=PhotoImage(file="wykr.gif")
   cv.create_image(0,0,anchor=NW,image=foto)
  except:
    pass


def cl():
  plt.clf()
  pole.delete(0,END)
  cv.delete('all')

def exit():
  plt.close()
  okno.destroy()

okno = Tk()
import tkinter.font
calcfont=tkinter.font.Font(font=("Courier", 10, "bold"))
ebg='#0099ff'
gadbg="#00ff00"
pole = Entry(okno)
pole.config(width=60, fg="white", bg=ebg, font=calcfont)
pole.grid(row=1, column=0, columnspan=8, pady=10)

button = Button(okno, text='1', command=onclick1, bg=gadbg, font=calcfont)
button.grid(row=3, column=0, sticky=EW)
button = Button(okno, text='2', command=onclick2, bg=gadbg, font=calcfont)
button.grid(row=3, column=1, sticky=EW)
button = Button(okno, text='3', command=onclick3, bg=gadbg, font=calcfont)
button.grid(row=3, column=2, sticky=EW)

button = Button(okno, text='4', command=onclick4, bg=gadbg, font=calcfont)
button.grid(row=4, column=0, sticky=EW)
button = Button(okno, text='5', command=onclick5, bg=gadbg, font=calcfont)
button.grid(row=4, column=1, sticky=EW)
button = Button(okno, text='6', command=onclick6, bg=gadbg, font=calcfont)
button.grid(row=4, column=2, sticky=EW)

button = Button(okno, text='7', command=onclick7, bg=gadbg, font=calcfont)
button.grid(row=5, column=0, sticky=EW)
button = Button(okno, text='8', command=onclick8, bg=gadbg, font=calcfont)
button.grid(row=5, column=1, sticky=EW)
button = Button(okno, text='9', command=onclick9, bg=gadbg, font=calcfont)
button.grid(row=5, column=2, sticky=EW)

button = Button(okno, text='0', command=onclick0, bg=gadbg, font=calcfont)
button.grid(row=6, column=0, sticky=EW)
button = Button(okno, text=',', command=onclickpoint, bg=gadbg, font=calcfont)
button.grid(row=6, column=1, sticky=EW)
button = Button(okno, text='=', command=onclickrow, bg=gadbg, font=calcfont)
button.grid(row=6, column=2, sticky=EW)

button = Button(okno, text='+', command=onclickpl, bg=gadbg, font=calcfont)
button.grid(row=2, column=4, sticky=EW)
button = Button(okno, text='-', command=onclickmin, bg=gadbg, font=calcfont)
button.grid(row=3, column=4, sticky=EW)
button = Button(okno, text='*', command=onclicktimes, bg=gadbg, font=calcfont)
button.grid(row=4, column=4, sticky=EW)
button = Button(okno, text='/', command=onclickdiv, bg=gadbg, font=calcfont)
button.grid(row=5, column=4, sticky=EW)
button = Button(okno, text='^', command=onclickpow, bg=gadbg, font=calcfont)
button.grid(row=6, column=4, sticky=EW)

"""
w = Label(okno, text="label")
w.grid(row=5, column=0)
"""


button = Button(okno, text='Rysuj', command=onclickRysuj, bg=gadbg, font=calcfont)
button.grid(row=2, column=5, sticky=EW)
button = Button(okno, text='Całka', command=onclickCalka, bg=gadbg, font=calcfont)
button.grid(row=3, column=5, sticky=EW)
button = Button(okno, text='M. Z.', command=onclickZero, bg=gadbg, font=calcfont)
button.grid(row=4, column=5, sticky=EW)
button = Button(okno, text='LOAD', command=onclickLoad, bg=gadbg, font=calcfont)
button.grid(row=5, column=5, sticky=EW)


button = Button(okno, text='EXIT', command=exit, bg='red', font=calcfont)
button.grid(row=6, column=5, sticky=EW)


button = Button(okno, text='C', command=cl, bg='red', font=calcfont)
button.grid(row=2, column=0, sticky=EW)
button = Button(okno, text='D', command=onclickdel, bg='#ff9933', font=calcfont)
button.grid(row=2, column=1, sticky=EW)
button = Button(okno, text='%', command=onclickmod, bg=gadbg, font=calcfont)
button.grid(row=2, column=2, sticky=EW)




cv = Canvas(okno, width=400, height=300)
cv["background"]="white"
cv["borderwidth"]=0
cv.config()
cv.grid(row=10, column=0, columnspan=9,pady=10)



okno.title('Kalkulator')
okno.mainloop()
