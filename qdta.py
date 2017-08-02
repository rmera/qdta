#!/usr/bin/env python

#A few functions to allow quick and interactive analysis of MD trajectories in PyMOL.
#I indent with tabs. Just use search/replace to change them to spaces if you prefer it that way.

import numpy as np
import matplotlib
from matplotlib.backends import backend_tkagg

#Although there is sitll TCL-threads mishandling (I get an error) this
#makes the script behave properly for me. I know these things are non-deterministic, so sorry if I make your
#nuclear facility explode.
#The code is from Dinar Abdullin on StackOverflow: https://stackoverflow.com/questions/29234531/create-a-new-tk-window-thread-in-the-pymol-session-and-output-matplotlib-graph-o
def _new_figure_manager(num, *args, **kwargs):
    # import pymol
    if pymol._ext_gui is None:
        return new_figure_manager(num, *args, **kwargs)
    backend_tkagg.show._needmain = False
    import tkinter as Tk
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, FigureManagerTkAgg
    FigureClass = kwargs.pop('FigureClass', Figure)
    print(kwargs)
    figure = FigureClass(*args, **kwargs)
    window = Tk.Toplevel(master=pymol._ext_gui.root)
    canvas = FigureCanvasTkAgg(figure, master=window)
    figManager = FigureManagerTkAgg(canvas, num, window)
    if matplotlib.is_interactive():
        figManager.show()
    return figManager

new_figure_manager = backend_tkagg.new_figure_manager
backend_tkagg.new_figure_manager = _new_figure_manager

import matplotlib.pyplot as plt
from chempy.models import Indexed
from chempy import Atom
from copy import deepcopy, copy


def unit_vec(vector):
    return vector / np.linalg.norm(vector)

#Some functions to be used with traj_plot
#Feel free to add functions to this. You can also just call plot_traj with whatever external function you want, as long as it takes N selections in a given state and returns one number.
class MD:
    #Returns the distance between selections. If one or both have more than one atom, the first one is considered for each of them.
    def Distance(sels,extra):
        if len(sels)<2:
            raise("Distance needs two selections!") #maybe this is too much?
        c1=np.array(sels[0].atom[0].coord)
        c2=np.array(sels[1].atom[0].coord)
        return np.linalg.norm(c2-c1)
    #Returns the angle  sels[0]---sels[1]---sels[2] in degrees. If the selections have more than one atom, the first one is taken for each of them. NOT TESTED!
    def Angle(sels,extra):
        if len(sels)<3:
            raise("Angle needs three selections!") 
        c1=np.array(sels[0].atom[0].coord)
        c2=np.array(sels[1].atom[0].coord)
        c3=np.arrray(sels[2].atom[0].coord)
        c1u=unit_vec(c1-c2)
        c3u=unit_vec(c3-c2)
        return (1 / 0.0174533)*np.arccos(np.clip(np.dot(c1u, c3u), -1.0, 1.0)) # this is in degrees
    #RMSD returns, well, the RMSD between two selections. They have to have the same amount of atoms, and
    #be superimposed already. The point of this is to use it in a trajectory, which you should have already superimposed with a reference
    #structure, which would be your first selection.
    #If you use this function with plot_traj you need to call it with reference_selections=1
    def RMSD(sels,extra):
        if len(sels)<2:
            raise("RMSD needs two selections!") 
        c=[]
        #There may very well be a more efficient way to do this.
        for i in range(len(sels[0].atom)):
            c.append(np.sum(np.array(sels[1].atom[i].coord)-np.array(sels[0].atom[i].coord))**2)
        return np.sqrt(np.sum(np.array(c))/len(c))


#A little helper. Returns an atom with the coordinates for the center of mass
#of the selection. Useful for when you want to use traj_plot on the center of mass of something, as opposed to on a single atom.
def com_object(selection):
    states=cmd.count_states(selection)
    outstates=[]
    for i in range(states):
        coords=cmd.centerofmass(selection,i+1)
        outstates.append(Indexed())
        outstates[-1].atom.append(Atom())
        outstates[-1].atom[0].name="GA"
        outstates[-1].atom[0].symbol="Ga"
        outstates[-1].atom[0].index=0
        outstates[-1].atom[0].resid=0
        outstates[-1].atom[0].resn=""
        outstates[-1].atom[0].coord=coords
    #    print("QL",outstates[-1].atom[0].coord)
     #   if i>=1:
     #       print("QL2",outstates[-2].atom[0].coord)
    for i in range(len(outstates)):
        cmd.load_model(outstates[i],"COM",1+i)
       # print(i, outstates[i], outstates[i].atom[0].coord)


def del_first(selection):
    states=cmd.count_states(selection)
    name=selection+"_minusfirst"
    for i in range(2,states+1):
        cmd.load_model(cmd.get_model(selection,i),name,i-1,discrete=0)

#takes a single state from the selection (the first state, by default) and returns it 
#as a new object. Mostly to be use as a reference for the MD.RMSD function.
def get_ref(selection,state=1):
    cmd.load_model(cmd.get_model(selection,state),"REF")

    

#traj_plot applies function to selections in all their states and plots the resulting values
#agains the state number.
#functions has to take a list of chempy objects (one per selection) for a given state plus a list of extra parameters and return ONE numerical value. 
#If reference_selections>0 first N selections will be given always in the first state
#the functions included in the MD class are ready to be used here, but you can easily define your own.
#The optional extra parameter is a list of additional parameters for the function called. It is empty by default.
#The optional function_info will be used as a label for the y-axis of theplot.
#Often the first state is just a reference, which is difficult to correctly take away in pymol (see our del_first for an attempt).
def traj_plot(selections,function, extra=[] ,reference_selections=-1,function_info="Function",skip_first=False):
    if selections=="":
        print('Usage:  traj_plot(selections,function, extra=[] ,reference_selections=-1,function_info="Function",skip_first=False)')
        exit(1) 
    #We kinda assume that all selections given have the same number of states!
    selections=selections.split(",") 
    # If reference_selections>0 the first reference_selections selections will only have one state
    #the last one should always have the full states.
    states=cmd.count_states(selections[-1])
    #print(states) ########################################################
    #it may be that states are counted from one, so one 
    y=[]
    x=[]
    ylabel=""
    for i in range(1,states+1):
        if i==1 and skip_first:
            continue
        objs=[]
        for j,v in enumerate(selections):
            k=i
            if j <= reference_selections-1:
                k=1 #we overwrite the state index for the first N=reference_selections selections
            objs.append(cmd.get_model(v,k))
        x.append(i)
        y.append(function(objs,extra))
    X=np.array(x)
    Y=np.array(y)
    #print("X and Y", X, Y)####################
    fig, ax = plt.subplots()
    ax.plot(X,Y,"b-") #blue line
    plt.xlabel("State")
    plt.ylabel(function_info)
    plt.show()





