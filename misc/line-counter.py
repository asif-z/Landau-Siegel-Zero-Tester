def count_lines(filename):
    try:
        with open(filename, 'r') as file:
            line_count = sum(1 for line in file)
        return line_count
    except FileNotFoundError:
        print(f"File '{filename}' not found.")
        return None

# Prompt the user for the file name
filename = input("Enter the file name: ")
line_count = count_lines(filename)

if line_count is not None:
    print(f"Number of lines in '{filename}': {line_count}")
