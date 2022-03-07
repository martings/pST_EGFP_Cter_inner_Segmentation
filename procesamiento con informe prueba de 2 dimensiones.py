#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 18:47:57 2019

@author: maxpower
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import skimage as ski
import pandas as pd
from skimage import data, img_as_float, measure, filters, io, morphology, segmentation, feature
from skimage.segmentation import (morphological_chan_vese,
                                  morphological_geodesic_active_contour,
                                  inverse_gaussian_gradient,
                                  checkerboard_level_set)


#%%

# Abrir la imagen
name="MAX_short2.tif"

print(name)
filenameG = '%s_EGFP.tif' %name
filenameR = '%s_RFP.tif' %name
imgR = io.imread(filenameR)
imgG = io.imread(filenameG)

fig, (ax0, ax1) = plt.subplots(nrows=1,
                                    ncols=2,
                                    figsize=(10, 4),
                                    sharex=True,
                                    sharey=True)

ax0.imshow(imgG)
ax0.set_title('EGFP')
ax0.axis('off')

ax1.imshow(imgR)
ax1.set_title('RFP')
ax1.axis('off')

#%%

def store_evolution_in(lst):
    """Returns a callback function to store the evolution of the level sets in
    the given list.
    """

    def _store(x):
        lst.append(np.copy(x))

    return _store


# Morphological ACWE
image = img_as_float(imgR)

# Initial level set
init_ls = checkerboard_level_set(image.shape, 10)
# List with intermediate results for plotting the evolution
evolution = []
callback = store_evolution_in(evolution)
ls = morphological_chan_vese(image, 50, init_level_set=init_ls, smoothing=5,
                             iter_callback=callback)


fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.2),
                                    sharex=True, sharey=True)




ax1.imshow(image, cmap="gray")
ax1.set_axis_off()
ax1.contour(ls, [0.5], colors='r')
ax1.set_title("Morphological ACWE segmentation", fontsize=12)

ax2.imshow(ls, cmap="gray")
ax2.set_axis_off()
contour = ax2.contour(evolution[10], [0.5], colors='g')
contour.collections[0].set_label("Iteration 10")
contour = ax2.contour(evolution[25], [0.5], colors='y')
contour.collections[0].set_label("Iteration 25")
contour = ax2.contour(evolution[50], [0.5], colors='r')
contour.collections[0].set_label("Iteration 50")
ax2.legend(loc="upper right")
title = "Morphological ACWE evolution"
ax2.set_title(title, fontsize=12)

fig.tight_layout()
plt.show()
#%%
plt.imshow(ls)
#%%

# nomenclar cada celula
msk = segmentation.clear_border(ls)
plt.imshow(msk);
plt.show(msk)

#%%
cells = measure.label(msk)
plt.imshow(cells);
# filtrar
for region in measure.regionprops(cells):
    
    if region.area < 500:
        cells[cells == region.label] = 0
        
cells, _, _ = segmentation.relabel_sequential(cells)
plt.imshow(cells)

#%%

plt.imshow(cells)
print(os.getcwd())

#%%
# Buscar nucleos (esto es trivial con un stain nuclear)
nuclei = np.zeros_like(cells)

for label in np.unique(cells):
    
    if label == 0:
        continue
    else:    
        cell_msk = cells==label
        threshold = filters.threshold_otsu(imgG)
        bw = morphology.closing(imgG > threshold)
        nuclei[np.logical_and(cell_msk, bw)] = label


plt.imshow(nuclei)

#%%
out = []
for label in np.unique(cells):
    
    if label == 0:
        continue
    else:    
        cell_msk = cells == label
        nuc_msk = nuclei == label
        out.append((label, np.sum(cell_msk), np.sum(imgR[cell_msk]), np.sum(nuc_msk), np.sum(imgG[nuc_msk])))

out = pd.DataFrame(out, columns=('label', 'cell_area', 'cell_integrated_intensity', 'nuclear_area', 'nuclear_integrated_intensity'))


#%%
plt.imshow(cells + nuclei)
out.to_csv('%s.csv' %name) 

#%%
fig, (ax0, ax1, ax2) = plt.subplots(nrows=1,
                                    ncols=3,
                                    figsize=(20, 10),
                                    sharex=True,
                                    sharey=True)

ax0.imshow(imgR)
ax0.set_title('original RFP')
ax0.axis('off')


ax1.imshow(imgG)
ax1.set_title('original GFP')
ax1.axis('off')


ax2.imshow(cells + nuclei)
ax2.set_title('procesada')
ax2.axis('off')


fig.savefig('%s.png' % name) 

#%%
from skimage.measure import label

# labeled contains one integer for each pixel in the image,
# where that image indicates the segment to which the pixel belongs
labeled = label(msk)
#%%
from skimage.measure import regionprops

# create array in which to store cropped articles
cropped_images = []
cropped_images_nuclei = []
cropped_images_original = []


# define amount of padding to add to cropped image
pad = 10


# for each segment number, find the area of the given segment.
# If that area is sufficiently large, crop out the identified segment.
for region_index, region in enumerate(regionprops(labeled)):
    if region.area < 1000:
        continue
# draw a rectangle around the segmented articles
# bbox describes: min_row, min_col, max_row, max_col
    minr, minc, maxr, maxc = region.bbox
# use those bounding box coordinates to crop the image
    cropped_images.append(cells[minr-pad:maxr+pad, minc-pad:maxc+pad])
    cropped_images_nuclei.append(imgG[minr-pad:maxr+pad, minc-pad:maxc+pad])
    cropped_images_original.append(imgR[minr-pad:maxr+pad, minc-pad:maxc+pad])
    
    #%%

import matplotlib.pyplot as plt
plt.imshow(cells+imgG)
#%%

for i in range(len(cropped_images)):
    plt.imshow(cropped_images[i]+cropped_images_nuclei[i])
    plt.figure(i+1)
plt.show()
#%%
os.chdir("/media/maxpower/Datos/Martin/Laboratorio/Microscopio/2019-06-17/citomCherry/short/stack/buenas/procesadas")
#%%
for i in range(len(cropped_images)):
    celu=cropped_images[i]
    nucleo=cropped_images_nuclei[i]
    original=cropped_images_original[i]
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(8, 2.5))
    yslice = celu.shape[0]//2
    if yslice == 0:
        continue
    #yslice = 50


    ax0.imshow(celu, vmin=image.min(), vmax=image.max())
    ax0.axhline(yslice, color='r', alpha=0.4)
    ax0.set_title('mascara mcherry')
    ax0.axis('off')

    ax1.imshow(nucleo)
    ax1.axhline(yslice, color='r', alpha=0.4)
    ax1.set_title('Itgb1 EGFP')
    ax1.axis('off')

    ax2.plot(celu[yslice], '0.5', label='mask')
    ax2.plot(nucleo[yslice], 'k', label='EGFP')
    #ax2.plot(original[yslice], 'r', label='mcherry')
    ax2.set_title('image slice')
    ax2.set_xticks([])
    ax2.legend()

    fig.tight_layout()
    plt.savefig('%s' %i + '%s' %name +  'horizontal.png')
plt.show()

#%%
for i in range(len(cropped_images)):
    celu=cropped_images[i]
    nucleo=cropped_images_nuclei[i]
    original=cropped_images_original[i]
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(8, 2.5))
    xslice = celu.shape[1]//2
    if xslice == 0:
        continue
    #yslice = 50
    ax0.imshow(celu, vmin=image.min(), vmax=image.max())
    ax0.axvline(xslice, color='r', alpha=0.4)
    ax0.set_title('mascara mcherry')
    ax0.axis('off')

    ax1.imshow(nucleo)
    ax1.axvline(xslice, color='r', alpha=0.4)
    ax1.set_title('Itgb1 EGFP')
    ax1.axis('off')

    ax2.plot(celu[:,xslice], '0.5', label='mask')
    ax2.plot(nucleo[:,xslice], 'k', label='EGFP')
    #ax2.plot(original[xslice], 'r', label='mcherry')
    ax2.set_title('image slice')
    ax2.set_xticks([])
    ax2.legend()

    fig.tight_layout()
    plt.savefig('%s' %i +'%s' %name + 'vertical.png')
plt.show()
#%%
print(nucleo[yslice])
#%%

#%%
for i in range(len(cropped_images)):
    celu=cropped_images[i]
    nucleo=cropped_images_nuclei[i]
    xslice = celu.shape[1]//2
    yslice = celu.shape[0]//2
    if xslice == 0:
        continue
    mcelu=np.array(celu[yslice])
    nuclesli=np.array(nucleo[yslice])
    asd=np.column_stack((mcelu,nuclesli))
    asd=pd.DataFrame(asd, columns=('cell_msk', 'EGFP' ))
    asd.to_csv('%s' %i +'%s' %name + 'horizontal.csv') 
print(asd)
#%%
for i in range(len(cropped_images)):
    celu=cropped_images[i]
    nucleo=cropped_images_nuclei[i]
    xslice = celu.shape[1]//2
    yslice = celu.shape[0]//2
    if xslice == 0:
        continue
    mcelu=np.array(celu[:,xslice])
    nuclesli=np.array(nucleo[:,xslice])
    asd=np.column_stack((mcelu,nuclesli))
    asd=pd.DataFrame(asd, columns=('cell_msk', 'EGFP' ))
    asd.to_csv('%s' %i +'%s' %name + 'vertical.csv') 
print(asd)
#%%