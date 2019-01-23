import glob
import sys

path = sys.argv[1]
pattern = path + "*.test"
files = glob.glob(pattern)
empty_files = 0

for filename in files:
  with open(filename, "r+") as f:
    contents = f.readlines()
    instance_name = filename.split("/")[-1].split(".")[0]
    if contents == []:
      empty_files += 1
      print instance_name, "empty"
      continue
    line0 = contents[0]
    if line0.split()[0] == 'inf':
      f.seek(0)
      f.write(instance_name + " inf inf inf inf inf inf\n")
      f.truncate()

#print "Empty files: ", empty_files
