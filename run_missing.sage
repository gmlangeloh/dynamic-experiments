import os
import subprocess

with open("missing.out", "r") as f:
  missing = f.readlines()

to_run = {}
for line in missing:
  words = line.split("/")
  algo_and_heuristic = words[0]
  algorithm = '-'.join(algo_and_heuristic.split('-')[:-1])
  heuristic = algo_and_heuristic.split('-')[-1]
  instance = "./instances/" + words[1].split(".")[0] + ".ideal"
  if (algorithm, heuristic) in to_run:
    to_run[(algorithm, heuristic)].append(instance)
  else:
    to_run[(algorithm, heuristic)] = [instance]

algorithms = {
    'static' : 'static',
    'cp' : 'caboara-perry',
    'random' : 'random',
    'caboara' : 'caboara',
    'perturb' : 'perturbation',
    'gs' : 'gritzmann-sturmfels',
    'simplex' : 'simplex',
    'regrets' : 'regrets',
    'gs-then-cp' : 'gs-then-cp',
    'simplex-then-cp' : 'simplex-then-cp'
}

def run(algorithm, heuristic):
  dirname = "./inst-" + algorithm + "-" + heuristic
  os.mkdir(dirname)
  for i in to_run[(algorithm, heuristic)]:
    subprocess.call('cp ' + i + ' ' + dirname, shell=True)
  args = [ algorithms[algorithm], heuristic, dirname ]
  subprocess.call(['./run_experiment.sh'] + args)
  #TODO Have to put results in correct directory
  target = algorithm + '-' + heuristic + '/'
  subprocess.call('mv *.test ' + target, shell=True)
  os.rmdir(dirname)

for (algorithm, heuristic) in to_run:
  run(algorithm, heuristic)
