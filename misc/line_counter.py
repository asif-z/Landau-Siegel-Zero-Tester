# Simple script that counts number of lines

fname = input("Enter the file name: ")

f = open(fname, "r")

counter = 0;

for line in f:
    counter += 1

f.close()

print("Line count: " + str(counter))
