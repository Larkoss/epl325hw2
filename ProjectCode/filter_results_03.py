def main():
    filename = "hw2_o3_hpc.txt"
    infile = open(filename, "r")
    for line in infile:
        elements = line.split()
        if len(elements) == 2:
            print()
            print(str(elements[1]) + "\t", end='')
        elif len(elements) == 3:
            print(elements[2] + "\t", end='')
        elif len(elements) == 0:
            i = 0
        else:
            i = 0
            #print()

    infile.close()


if __name__ == '__main__':
    main()
