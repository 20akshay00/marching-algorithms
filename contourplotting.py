import matplotlib.pyplot as plt

def contourplot():
    file1 = open("edges.txt", 'r')
    lines = file1.readlines()

    for line in lines:
        line = line.strip()
        line = line.split()
        line = [float(i) for i in line]

        plt.plot([line[0], line[2]], [line[1], line[3]], 'bo-')
        print(line)
    plt.grid()
    plt.show(block=True)

    file1.close()
    
contourplot()

    


