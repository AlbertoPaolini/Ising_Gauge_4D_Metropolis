import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def write_to_file(filename, data):
    with open(filename, 'a') as f:
        for line in data:
            f.write(line + "\n")

write_to_file("Data.txt", ["Block_Size Mean Standard_Deviation Variance Beta Standard_Deviation_of_Variance Variance_Unblocked"])  # Formato CSV

print("Insert the number of dataset:")
n1 = int(input())

print("How many data-points do you want to reject because of the thermalization?")
reject = int(input())

print("How many data-points do you want to enter in each block?")
max_size = int(input())
folder_path = "Energy_data"

file_names = [f"energy_data{i}.txt" for i in range(1, n1+1)]

plt.figure()

for file_name in file_names:
    file_path = os.path.join(folder_path, file_name)
    
    if os.path.isfile(file_path):  
        print(f"Processing {file_name}...")
        
        data = np.loadtxt(file_path, usecols=(1, 2))
        data_values = data[:, 0]
        beta = data[0, 1]
        
        data_values = data_values[reject:]

        variance = []
        std = []
        average = []
        n = []
        variance_unblocked = np.var(data_values)
        
        max_valid_size = min(max_size, len(data_values))  

        for size in range(1, max_valid_size + 1, 1):
            
            blocks_average = []
            variance_average = []
            num_blocks = len(data_values) // size
            
            for i in range(num_blocks):
                block = data_values[i * size : (i + 1) * size]
                blocks_average.append(np.mean(block))
                variance_average.append(np.var(block))
            
            if blocks_average:
            
                block_mean = np.mean(blocks_average)
                block_variance = np.var(blocks_average)
                block_std = np.sqrt( block_variance / num_blocks )
                
                variance_mean = np.mean( variance_average )
                variance_variance = np.var( variance_average )
                variance_std = np.sqrt( variance_variance / num_blocks )
            
            else:
                block_mean = block_variance = block_std = 0
            
            variance.append(block_variance)
            average.append(block_mean)
            std.append(block_std)
            n.append(size)
            
        plt.plot(n, std, marker='o', label=f"Beta = {beta}")

        result_line = f"{size} {block_mean} {block_std} {block_variance} {beta:.5f} {variance_std} {variance_unblocked}"
        write_to_file("Data.txt", [result_line])

plt.xlabel('Data-points per block', fontsize=14)
plt.ylabel('Standard deviation', fontsize=14)
formatter = FuncFormatter(lambda x, _: f'{x:.0e}')  # Notazione scientifica con due decimali
plt.gca().yaxis.set_major_formatter(formatter)

plt.title('Blocking method for the density of energy', fontsize=14)
plt.legend()

plt.show()
