#authors:Erik Ernesto Ocegueda Sambrano y Paola Enríquez Reyes
#date:6/10/2023

#INTRODUCTION
#NAME: FAKESKY
#PURPOSE:
#      MASKOUT THE STARS ON THE EXTRAORDINARY STARS.
#      THE PROCEDURE CAN BE MODIFY TO MASKOUT THE ORDINARY STARS.
#      THE MASKOUT DEPENDS ON THEIR BESTAPP VALUE.
#      MASK USES THE MEAN VALUE OF THE BACKGROUND.
#      TAKES MAXIMUM BESTAPP 10, THIS CAN BE MODIFY.
#
# EXPLANATION:
#        
#	
# CALLING SEQUENCE:
#     FAKESKY(image.fits, text_file)
#
# INPUTS:
#       imagen_fits - image of reference to produce the fake sky
#
#	    calcpol - file with positions in
#                  columns: ID   X       Y       BestApp	 
#
# OPTIONAL INPUTS:
#     
# OUTPUTS:
#      
# 
#  EXAMPLE:
#
#      FAKESKY,'hd110cv0001.fits','calcpol_new .out'


import matplotlib.pyplot as plt
import astropy.io.fits as fits #pip install astropy
import numpy as np

fig, ax = plt.subplots()

#Path of image FITS
ruta_imagen = 'bias_flat_HD181474_L0_0avgcomb_1.fits'

#Read image FITS
imagen_fits = fits.open(ruta_imagen)[0]
datos_imagen = imagen_fits.data
datos_hdr = imagen_fits.header

#Read calcpol
data = []
with open('calcpol.out', 'r') as file:
    next(file)
    for line in file:
        values = line.split()
        data.append([int(values[0]), float(values[1]), float(values[2]), int(values[3])])
if (len(data) == 0): #Stops with a warning if there is no data found in the calcpol
    print("WARNING -NO DATA FOUND in calcpol.out-")
else:
        
    # Group data by BestApp
    grouped_data = {}
    for row in data:
        app = row[3]
        if app in grouped_data:
            grouped_data[app].append(row)
        else:
            grouped_data[app] = [row]

    Xval = []
    Yval = []
    sizes = []
    for app, group in grouped_data.items():
        if app <= 10:  #BESTAPP MAXIMUM
            X = [row[1] for row in group]
            Y = [row[2] for row in group]
            size = [row[3] for row in group]
            Xval.append(X)
            Yval.append(Y)
            sizes.append(size)
    cont = 0
    mean_value = np.mean(datos_imagen)

    i = 0
    j = 0

    i = 0
    while i < len(Xval):
        j = 0
        while j < len(Xval[i]):
            size = sizes[i][j]
            if size <= 10:  #BESTAPP MAXIMUM
                if j % 2 != 0:
                    x_pos = int(Xval[i][j])
                    y_pos = int(Yval[i][j])
                    radius = size * 1.5  # Radius circle
                    for y_offset in range(-int(radius), int(radius) + 1):
                        for x_offset in range(-int(radius), int(radius) + 1):
                            x = x_pos + x_offset
                            y = y_pos + y_offset
                            #Checks if the pixel is inside the circle
                            if (
                                0 <= x < datos_imagen.shape[1]
                                and 0 <= y < datos_imagen.shape[0]
                                and (x - x_pos) ** 2 + (y - y_pos) ** 2 <= radius ** 2
                            ):
                                datos_imagen[y][x] = mean_value
                j += 1
            else:
                half_size = size // 2
                for _ in range(size):
                    if j % 2 != 0:
                        x_pos = int(Xval[i][j])
                        y_pos = int(Yval[i][j])
                        radius = half_size * 1.5  # Radius circle
                        for y_offset in range(-int(radius), int(radius) + 1):
                            for x_offset in range(-int(radius), int(radius) + 1):
                                x = x_pos + x_offset
                                y = y_pos + y_offset
                                #Checks if the pixel is inside the circle
                                if (
                                    0 <= x < datos_imagen.shape[1]
                                    and 0 <= y < datos_imagen.shape[0]
                                    and (x - x_pos) ** 2 + (y - y_pos) ** 2 <= radius ** 2
                                ):
                                    datos_imagen[y][x] = mean_value
                    j += 1
        i += 1



    ax.plot()
    ax.imshow(datos_imagen, vmin = np.min(datos_imagen), vmax = 255, origin = 'lower')
    plt.title('Image')

    plt.show()


    # Saves file FITS in your PC
    hdu = fits.PrimaryHDU(datos_imagen, header = datos_hdr)
    hdu.writeto('fakesky.fits', overwrite=True)

    #Final
    print('FAKESKY run well: output fakesky.fits')
