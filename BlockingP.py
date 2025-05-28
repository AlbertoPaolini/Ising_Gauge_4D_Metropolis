import os
import numpy as np
import matplotlib.pyplot as plt

def write_to_file(filename, data):
    with open(filename, 'a') as f:
        for line in data:
            f.write(line + "\n")

write_to_file("Polyakov.txt", ["Block_Size Mean Standard_Deviation Variance Beta Standard_Deviation_of_Variance Variance_Unblocked"])  

print("How many dataset do you have?")
n1 = int(input())

print("How many data-points do you want to enter in each block?")
max_size = int(input())
folder_path = "Polyakov_data"

file_names = [f"polyakov_data{i}.txt" for i in range(1, n1+1)]

plt.figure()

for file_name in file_names:
    file_path = os.path.join(folder_path, file_name)
    
    if os.path.isfile(file_path):  
        print(f"Processing {file_name}...")
        
        data = np.loadtxt(file_path, usecols=(0, 1))
        data_values = data[:, 0]
        beta = data[0, 1]

        variance_unblocked = np.var(data_values, ddof=1)
        print(variance_unblocked)
        print(np.mean(data_values**2)-np.mean(data_values)**2)

        variance = []
        std = []
        average = []
        n = []
        
        max_valid_size = min(max_size, len(data_values))  

        for size in range(1, max_valid_size + 1, 2):
            
            block_average = []
            block_variance = []
            
            num_blocks = len(data_values) // size 

            for i in range(num_blocks):

                block = data_values[i * size : (i + 1) * size]
                block_average.append(np.mean(block))
                block_variance.append( np.var(block) )
                
            
            if block_average:
                blocks_mean = np.mean(block_average)
                blocks_variance = np.var(block_average)
                blocks_std = np.sqrt(blocks_variance / num_blocks)

                variance_variance = np.var(block_variance)
                variance_std = np.sqrt( variance_variance/ num_blocks)
            
            else:
                block_mean = block_variance = block_std = 0
            
            variance.append(blocks_variance)
            average.append(blocks_mean)
            std.append(blocks_std)
            n.append(size)
        plt.plot(n, std, marker='o', label=f"Beta = {beta}")

        result_line = f"{size} {blocks_mean} {blocks_std} {blocks_variance} {beta:.5f} {variance_std} {variance_unblocked}"
        write_to_file("Polyakov.txt", [result_line])

plt.xlabel('Data-points per block', fontsize=14)
plt.ylabel('Standard deviation', fontsize=14)
plt.title('Blocking method for the Polyakov loop', fontsize=16)
plt.legend()

plt.show()
