import os,sys
import numpy as np
import torch

#from normalization_constants import mean, std

def makeImage(prong_vv):
  plane0pix_row = []
  plane0pix_col = []
  plane0pix_val = []
  plane1pix_row = []
  plane1pix_col = []
  plane1pix_val = []
  plane2pix_row = []
  plane2pix_col = []
  plane2pix_val = []
  raw_plane0pix_row = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_col = np.zeros(prong_vv[3].size(), dtype=int)
  raw_plane0pix_val = np.zeros(prong_vv[3].size(), dtype=float)
  raw_plane1pix_row = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_col = np.zeros(prong_vv[4].size(), dtype=int)
  raw_plane1pix_val = np.zeros(prong_vv[4].size(), dtype=float)
  raw_plane2pix_row = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_col = np.zeros(prong_vv[5].size(), dtype=int)
  raw_plane2pix_val = np.zeros(prong_vv[5].size(), dtype=float)
  if prong_vv[0].size()>=10:
    for i, pix in enumerate(prong_vv[0]):
      if pix.inCrop:
        plane0pix_row.append(pix.row)
        plane0pix_col.append(pix.col)
        plane0pix_val.append(pix.val)
  if prong_vv[1].size()>=10:
    for i, pix in enumerate(prong_vv[1]):
      if pix.inCrop:
        plane1pix_row.append(pix.row)
        plane1pix_col.append(pix.col)
        plane1pix_val.append(pix.val)
  if prong_vv[2].size()>=10:
    for i, pix in enumerate(prong_vv[2]):
      if pix.inCrop:
        plane2pix_row.append(pix.row)
        plane2pix_col.append(pix.col)
        plane2pix_val.append(pix.val)
  if prong_vv[3].size()>=10:
    for i, pix in enumerate(prong_vv[3]):
      raw_plane0pix_row[i] = pix.row
      raw_plane0pix_col[i] = pix.col
      raw_plane0pix_val[i] = pix.val
  if prong_vv[4].size()>=10:
    for i, pix in enumerate(prong_vv[4]):
      raw_plane1pix_row[i] = pix.row
      raw_plane1pix_col[i] = pix.col
      raw_plane1pix_val[i] = pix.val
  if prong_vv[5].size()>=10:
    for i, pix in enumerate(prong_vv[5]):
      raw_plane2pix_row[i] = pix.row
      raw_plane2pix_col[i] = pix.col
      raw_plane2pix_val[i] = pix.val
  plane0pix_row = np.array(plane0pix_row, dtype=int)
  plane0pix_col = np.array(plane0pix_col, dtype=int)
  plane0pix_val = np.array(plane0pix_val, dtype=float)
  plane1pix_row = np.array(plane1pix_row, dtype=int)
  plane1pix_col = np.array(plane1pix_col, dtype=int)
  plane1pix_val = np.array(plane1pix_val, dtype=float)
  plane2pix_row = np.array(plane2pix_row, dtype=int)
  plane2pix_col = np.array(plane2pix_col, dtype=int)
  plane2pix_val = np.array(plane2pix_val, dtype=float)
  
  image = np.zeros((6,512,512))
  image[0, plane0pix_row, plane0pix_col] = plane0pix_val
  image[2, plane1pix_row, plane1pix_col] = plane1pix_val
  image[4, plane2pix_row, plane2pix_col] = plane2pix_val
  image[1, raw_plane0pix_row, raw_plane0pix_col] = raw_plane0pix_val
  image[3, raw_plane1pix_row, raw_plane1pix_col] = raw_plane1pix_val
  image[5, raw_plane2pix_row, raw_plane2pix_col] = raw_plane2pix_val
  #image_t = torch.from_numpy(image).float()
  #norm = transforms.Normalize(mean, std)
  #image_t = norm(image_t).reshape(1,6,512,512)
  #return torch.clamp(image_t, max=4.0),image
  return image
