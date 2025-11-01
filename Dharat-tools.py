import math
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal 
from scipy.interpolate import griddata
from scipy import interpolate
from scipy.interpolate import LinearNDInterpolator
import warnings

# Ignore all warnings
warnings.filterwarnings("ignore")

print('This tool is for landslide slip surface estimation using spline interpolation method. For more informaiton visit: https://github.com/jsidhu45/Dharat-tools')

time.sleep(3)

data_dir = './input_data/'
output_dir = './output_files/'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
results_dir = './results/'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

dem_path = input("Enter the name of DEM (.tif) file: ")
full_dem_path = os.path.join(data_dir,dem_path)

res=float(input("Input the DEM resolution (also resample the DEM): "))

boundary_path = input("Enter the name of landslide boundary (.shp) file: ")
full_boundary_path = os.path.join(data_dir, boundary_path)

it=int(input("Total number of iterations: "))

result_output=float(input("Enter number of iterations after which it save the interpolated excel file: "))

dem_org = gdal.Open(full_dem_path)

dem_org=gdal.Warp(output_dir+'dem_resampled.tif', dem_org, xRes=res, yRes=res, resampleAlg="bilinear")

dem_band=dem_org.GetRasterBand(1)
dem_array=dem_band.ReadAsArray()
gt=dem_org.GetGeoTransform()

r,c=np.shape(dem_array)
dem_resamp=dem_org


box_xyz=gdal.Translate(output_dir+'box_xyz.xyz', dem_org)

df_box=pd.read_csv(output_dir+'box_xyz.xyz', sep=" ", header=None)
df_box.columns=['x', 'y', 'z']
df_box1=df_box
df_box_real=df_box

#clip dem of the landslide area using shp file
 
dem_clip=gdal.Warp(output_dir+'DEM_clip.tif',dem_org, cutlineDSName= full_boundary_path, dstNodata=np.nan)
clip_arr=dem_clip.GetRasterBand(1).ReadAsArray()


ls_inn_z=[]
for i in clip_arr.flatten():
    if i>0:
        ls_inn_z.append(i)


#finidng point data lying inside the landslide shape file

inn_xyz=gdal.Translate(output_dir+'inn_xyz.xyz', dem_clip)
df_inn_xyz=pd.read_csv(output_dir+'inn_xyz.xyz', sep=" ", header=None)
df_inn_xyz.columns=['x', 'y', 'z']

ls_inn=[]
for i in df_inn_xyz.to_numpy().tolist():
    if i[2] >0:
        ls_inn.append(i)
        
df_inn=pd.DataFrame(ls_inn)
df_inn.columns=['x', 'y', 'z']

df_inn_real=df_inn


#counting number of points in x and y direction lying outside the landslide shape file 

count_x=[]
count_x1=[]
for i in df_box_real['x'].to_numpy().tolist():
    if i < np.min(df_inn_real['x']):
        count_x.append(i)
        count_x1.append(i)
    elif i > np.max(df_inn_real['x']):
        count_x.append(i)

sub_c=len(np.unique(count_x))


count_y=[]
count_y1=[]
for i in df_box_real['y'].to_numpy().tolist():
    if i < np.min(df_inn_real['y']):
        count_y.append(i)
    elif i > np.max(df_inn_real['y']):
        count_y.append(i)
        count_y1.append(i)
        
sub_r=len(np.unique(count_y))



#Spline SLBL calculations
df_inn=round(df_inn,3)                                                  #can be replaced
df_box=round(df_box,3)                                                  #can be replaced

df_box=df_box.sort_values(by = ['x', 'y'], ascending = [True, False])

df_inn=df_inn.sort_values(by = ['x', 'y'], ascending = [True, False])

df_box_arr = df_box.to_numpy()


df_box_l = [i.tolist() for i in df_box_arr]

df_inn_l = df_inn.to_numpy().tolist()

zbox_l=df_box['z'].to_numpy().tolist()
zbox_arr=np.array(zbox_l)


mat_boxz=zbox_arr.reshape(c,r).T


xbox_l=df_box['x'].to_numpy().tolist()
xbox_arr=np.array(xbox_l)
mat_boxx=xbox_arr.reshape(c,r).T

ybox_l=df_box['y'].to_numpy().tolist()
ybox_arr=np.array(ybox_l)
mat_boxy=ybox_arr.reshape(c,r).T


mat_boxz=zbox_arr.reshape(c,r).T

mat_boxz1=np.zeros((r,c))

for i in range(c):
    for j in range(r):
        mat_boxz1[j,i]= mat_boxz[j,i]

avg1g=[]
avg1=[]
xyz1=[]



start=time.time()

itr=0
for t in range(it):
    for j in range(r-sub_r):
        for i in range(c-sub_c):
            indi=len(np.unique(count_x1))
            indj=len(np.unique(count_y1))
            z1=mat_boxz[j+indj,i+indi]
            x1=mat_boxx[j+indj,i+indi]
            y1=mat_boxy[j+indj,i+indi]
            xyz=[x1,y1,z1]
            
            for d in df_inn_l:
                if d == xyz:
                    rr_z = mat_boxz1[j+indj].tolist()
                    rr_x = mat_boxx[j+indj].tolist()
                    del rr_z[i+indi]
                    del rr_x[i+indi]
                    
    
                    f=interpolate.interp1d(rr_x,rr_z, 'cubic')                   #spline interpolation along x-axis
                    p = f(mat_boxx[j+indj, i+indi])
                    pt1 = p
                    
                    cc_z1 = mat_boxz1[:,i+indi].tolist()
                    cc_y1 = mat_boxy[:,i+indi].tolist()
                    del cc_z1[j+indj]
                    del cc_y1[j+indj]

                    f1=interpolate.interp1d(cc_y1,cc_z1, 'cubic')               #spline interpolation along y-axis
                    p1 = f1(mat_boxy[j+indj, i+indi])
                    pt2 = p1
                    
                    avg = (pt1+pt2)*0.5
                    if avg < z1: 
                        mat_boxz1[j+indj,i+indi] = avg
                    else:
                        mat_boxz1[j+indj,i+indi] = z1
    
    itr = itr+1
    if (t+1) % result_output == 0:
    	pd.DataFrame(mat_boxz1).to_excel(results_dir+ 'interpolated'+str(itr)+'.xlsx')
    end=time.time()
    
print("Time per iteration:", (end-start)/itr, "seconds")
pd.DataFrame(mat_boxz1).to_excel(results_dir+ 'interpolated.xlsx') 
pd.DataFrame(mat_boxz).to_excel(results_dir+ 'original.xlsx') 

print(f"The interpolated results files (.xlsx) are saved in {os.getcwd()}/results")

print('For more infomation read and cite : Singh, J and Sepulveda, SA., (2025) Assessing the landslide failure surface depth and volume: A new spline interpolation method. Engineering Geology. https://doi.org/10.1016/j.enggeo.2025.108319')
