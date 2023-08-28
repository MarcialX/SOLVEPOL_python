# Open the file in read mode
import glob
import os
import random
import math 

random.seed(0)

# function get_table reads the text file and outputs the desired table we wish to read data from. 
def get_table(file_path):
    with open(file_path, 'r') as file:
        table_started = False
        table_data = []

        # Iterate over each line in the file
        for line in file:
            # Check if the line starts with "PAIR"
            if line.startswith("PAIR"):
                column_names = line.strip().split()[1:]
                table_started = True
                continue

            # Check if the line is empty or contains only whitespace
            if not line.strip():
                table_started = False

            # If table has started, extract the data
            if table_started:
                data = line.strip().split()
                table_data.append(data)
        return table_data

# This part of the code reads the table obtained through above function, and stores values of fluxes
     # and errors in different lists

flux_ord = [] # stores all the ordinary flux values from all the files 
flux_ex = [] # stores all the extra_ordinary flux values from all the files 
error_ord = [] # stores all the ordinary error values from all the files 
error_ex = [] # stores all the extra_ordinary error values from all the files 
error_sq_e = [] # extra_ordinary error square values
error_sq_o = [] # ordinary error square values


# please change path to the folder in which .phot files are there. 
for file in sorted(glob.glob("/Users/Mansoor/pre_process/astro/*.phot")): # read .phot files 
    filename = os.path.basename(file) # obtain base name from path
    filename = os.path.splitext(filename)[0] # get filename from basename
    table = get_table(file) # function to read the table
    error_o = []
    error_e = []
    flux_o = []
    flux_e = []
    err_o = []
    err_e = []
    for count, i in enumerate(table): # iterate through the table
        if count % 2 == 1:
            flux_o.append([float(x) for x in i[8:18]])
            error_o.append(i[18:28])
            z = [float(x) for x in i[18:28]]
            err_o.append([x**2 for x in z])
        else:
            flux_e.append([float(x) for x in i[8:18]])
            error_e.append(i[18:28])
            z = [float(x) for x in i[18:28]]
            err_e.append([x**2 for x in z])
    flux_ord.append(flux_o)
    flux_ex.append(flux_e)
    error_ord.append(error_o)
    error_ex.append(error_e)
    error_sq_o.append(err_o)
    error_sq_e.append(err_e)
     
 
# This part of the code prepares summation of ordinary/extra_ordinary flux values, and error values. 
add_o = [] # stores summation of ordinary flux values
add_e = [] # stores summation of extra ordinary flux values
sum_err_o = [] # stores summation of ordinary error values
sum_err_e = [] # stores summation of extra ordinary error values
for j in range(len(flux_ord[0])):
    lists_odd = [flux_ord[0][j], flux_ord[1][j], flux_ord[2][j], flux_ord[3][j], flux_ord[4][j], flux_ord[5][j],
                 flux_ord[6][j], flux_ord[7][j]]
    lists_even = [flux_ex[0][j], flux_ex[1][j], flux_ex[2][j], flux_ex[3][j], flux_ex[4][j], flux_ex[5][j],
                  flux_ex[6][j], flux_ex[7][j]]
    error_odd = [error_sq_o[0][j], error_sq_o[1][j], error_sq_o[2][j], error_sq_o[3][j], error_sq_o[4][j], error_sq_o[5][j],
                 error_sq_o[6][j], error_sq_o[7][j]]
    error_even = [error_sq_e[0][j], error_sq_e[1][j], error_sq_e[2][j], error_sq_e[3][j], error_sq_e[4][j], error_sq_e[5][j],
                  error_sq_e[6][j], error_sq_e[7][j]]
     
    sum_o = [sum(x) for x in zip(*lists_odd)]
    sum_e = [sum(x) for x in zip(*lists_even)]
    error_odd = [math.sqrt(sum(x)) for x in zip(*error_odd)]
    error_even = [math.sqrt(sum(x)) for x in zip(*error_even)]
    add_o.append(sum_o)
    add_e.append(sum_e)
    sum_err_o.append(error_odd)
    sum_err_e.append(error_even)
     
# This part of the code computes all the parameters ansd stores in text file 
angles = [0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5] # define anlges 
myfile = open('result.txt', 'w') # define text file for the results 
myfile.write("PAIR 1" + '\n') # print on the top in text file 
print("len", len(flux_ord))
for j in range(len(flux_ord[0])): #  Here j runs for number equal to n umber of files. If files are 8 then j=8 
    for k in range(len(flux_ord[0][0])): # here k runs for 10 values since there are 10 flux values
        result = []
        dzz = []
        q = 0
        u = 0
        dq = 0
        du = 0
        for l in range(len(flux_ord)): #  here l=number of files if files are 8
            z = 0
            z1 = flux_ord[l][j][k] - flux_ex[l][j][k] * (add_o[j][k] / add_e[j][k])
            z2 = flux_ord[l][j][k] + flux_ex[l][j][k] * (add_o[j][k] / add_e[j][k])
            result.append(round(z1 / z2, 5))
            ak = float(add_e[j][k])/float(add_o[j][k])
            dak = math.sqrt(sum_err_e[j][k]**2 + ak**2 * sum_err_o[j][k]**2)/abs(float(add_o[j][k]) )
            dz = math.sqrt((flux_ord[l][j][k]**2 * ak**2 * float(error_ex[l][j][k])**2)+ 
                           (flux_ex[l][j][k]**2 * ak**2 * float(error_ord[l][j][k])**2)
                           +(flux_ex[l][j][k]**2 * flux_ord[l][j][k]**2 * dak**2))/ (flux_ex[l][j][k]+flux_ord[l][j][k]*ak)**2
            dzz.append(round(2*dz, 5))
            q+= 2* (round(z1 / z2, 5)*math.cos(4*angles[l] * math.pi / 180))/len(flux_ord)
            u+= 2* (round(z1 / z2, 5)*math.sin(4*angles[l] * math.pi / 180))/len(flux_ord)
            dq += 2 * math.sqrt(dz**2 * abs(math.cos(4*angles[l] * math.pi / 180))) / len(flux_ord)
            du += 2 * math.sqrt(dz**2 * abs(math.sin(4*angles[l] * math.pi / 180))) / len(flux_ord)
        p = math.sqrt(q**2 + u**2)
        dp = round(math.sqrt((q**2 * round(dq, 5)**2) + (u**2 *round(du, 5)**2)) / p, 5)
        theta = math.degrees(math.atan(u/q))
        z_sum = sum(i*i for i in result)
        sigma = math.sqrt(((z_sum* 0.25) - (q **2 + u**2))/ 6)
        dtheta = 28.65*sigma / p
       
        # storing values to the text file
        myfile.write("z values: Aperture {}".format(k + 1)+ (" " * 4) + str(result) + '\n')
        myfile.write("dz values: Aperture {}".format(k + 1)+ (" " * 4) + str(dzz) + '\n')
        myfile.write("q" + (" " * 13) + "dq" + (" " * 12) + "u" + (" " * 8) + "du" + (" " * 10) + "p" + (" " * 10) + "sigma" + (" " * 10) + "dp" + (" " * 10) + "theta" + (" " * 10) + "dtheta" + '\n')
        myfile.write(str(round(q, 5))+ (" " * 4) + str(round(dq, 5))+ (" " * 4) +  str(round(u, 5))+ (" " * 4) + str(round(du, 5))+ (" " * 4) + str(round(p, 5))+ (" " * 7) + str(round(sigma, 5)) + (" " * 7) + str(round(dp, 5)) + (" " * 6) + str(round(theta, 5)) + (" " * 7) + str(round(dtheta, 5))+'\n')
        myfile.write('\n')
    if j+2< 54:
        myfile.write('\n' + "PAIR {}".format(j + 2) + '\n')
         

myfile.close()
