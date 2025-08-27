import os
import csv
import re

def search_fail_in_rank_files(directory='.'):

    f = open("output.txt","w")

    keyword = "fail"
    keyword_lower = keyword.lower()

    print(f"Searching for '{keyword}' in output_rank_*.csv files under: {directory}")
    print("-" * 60)

    pattern = re.compile(r'^output_rank_\d+\.csv$')

    for root, _, files in os.walk(directory):
        for file in files:
            if pattern.match(file):
                print(file)
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, mode='r', encoding='utf-8') as csvfile:
                        reader = csv.reader(csvfile)
                        for row_num, row in enumerate(reader, start=1):
                            if any(keyword_lower in str(cell).lower() for cell in row):
                                # print(f"Match in: {file_path}")
                                # print(f"Row {row_num}: {row}")
                                f.write(row[0]+"\n")
                                # print("-" * 60)
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")

# Run the search
if __name__ == '__main__':
    search_dir = input("Enter directory to search (leave blank for current): ").strip()
    if not search_dir:
        search_dir = '.'
    search_fail_in_rank_files(search_dir)
