import sys
import os

read_path = input("Enter the path of your (paired-end) reads: ")
read_ext = ".gz"

assert os.path.exists(read_path), "Check path - unable to find directory: "+str(read_path)
read_list = sorted([i for i in os.listdir(read_path) if os.path.splitext(i)[1] == read_ext])
if len(read_list) > 0 :
	print("Reads found:")	
	for i in read_list :
		print(str(i))
	with open('sample_ids.txt', 'w') as f:
		f.write("sample_id\n")		
		for item in read_list:
			f.write("%s\n" % item)
else :
	print("No files found ending in \".fastq.gz\" in directory \""+str(read_path)+ "\"")
#stuff you do with the file goes here
#f.close()

database_path = input("Enter the path of your gene databases: ")
database_ext = ".prepareref"

assert os.path.exists(database_path), "Check path - unable to find directory: "+str(database_path)
database_list = sorted([i for i in os.listdir(database_path) if os.path.splitext(i)[1] == database_ext])
if len(database_list) > 0 :
	print("Databases found:")	
	for i in database_list :
		print(str(i))
	with open('gene_dbs.txt', 'w') as f:
		f.write("gene_database\n")		
		for item in database_list:
			f.write("%s\n" % item)
else :
	print("No files found ending in \".prepareref\" in directory \""+str(read_path)+ "\"")
#stuff you do with the file goes here
#f.close()

MLSTDB_path = input("Enter the path of your mlst databases: ")
MLSTDB_ext = ".mlstdb"

assert os.path.exists(database_path), "Check path - unable to find directory: "+str(database_path)
MLSTDB_list = sorted([i for i in os.listdir(MLSTDB_path) if os.path.splitext(i)[1] == MLSTDB_ext])
if len(MLSTDB_list) > 0 :
	print("MLST databases found:")	
	for i in MLSTDB_list :
		print(str(i))
	with open('mlst_dbs.txt', 'w') as f:
		f.write("mlst_database\n")		
		for item in MLSTDB_list:
			f.write("%s\n" % item)
else :
	print("No files found ending in \".mlstdb\" in directory \""+str(MLSTDB_PATH)+ "\"")
#	print("No files found ending in \".prepareref\" in directory \""+str(database_path)"\"")
#stuff you do with the file goes here
#f.close()
