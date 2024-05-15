from subprocess import check_output, CalledProcessError
import os

# Path to the model to be orchestrated
template_path = "moving_boundaries_grid3D_4scaling.py"
temp_model_path = "moving_boundaries_grid3D_4scaling_tmp.py"

# Load the file to be executed
with open(template_path, "r") as template_file:
  model_template = template_file.read()

# Validate template contains required symbols once each
if model_template.count("<N_CELLS>") != 1:
    raise Exception("Model template should contain a single occurrence of replacement symbol '<N_CELLS>'")
if model_template.count("<ECM_N>") != 1:
    raise Exception("Model template should contain a single occurrence of replacement symbol '<ECM_N>'")

# For each configuration
N_CELLS_sweep = [1000, 10000, 100000, 1000000]
ECM_N_sweep = range(10,101,10)
for N_CELLS in N_CELLS_sweep:
    for ECM_N in ECM_N_sweep:
        # Perform replacement of the two variables
        # @note These replacements are unstable if the template changes
        model_src = model_template.replace("<N_CELLS>", "%d"%(ECM_N), 1).replace("<ECM_N>", "%d"%(ECM_N), 1)
        with open(temp_model_path, "w") as model_file:        
            f.write(model_src)
        # Execute (might require python3)
        try:
            output = check_output(["python", temp_model_path], stderr=STDOUT).decode()
        except CalledProcessError as e:
            output = e.output.decode()
        # Collect/validate performance results
        

# Cleanup
if os.path.exists(temp_model_path)
    os.remove(temp_model_path)