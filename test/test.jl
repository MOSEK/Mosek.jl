using Mosek

env = makeenv()

task = maketask(env)

model_filename = "/home/idunning/mosek/7/tools/examples/data/25fv47.mps"
readdata(task, model_filename)

optimize(task)

output_filename = "/home/idunning/out.txt"
writedata(task, output_filename)
