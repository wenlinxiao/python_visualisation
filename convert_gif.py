#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Combining multiple .png figures into a .gif

Created on Wed Nov 11 03:47:58 2020

@author: Wenlin Xiao
"""

import imageio as im
import os

path="./gif/"
im_fns=sorted(os.listdir(path))
gif_imgs=[]
for imfn in im_fns:
    gif_imgs.append(im.imread(path+imfn))
im.mimsave("Rossby.gif",gif_imgs,fps=1)


