from subprocess import check_output, CalledProcessError
import os, re, csv
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

# Output regex captures
sim_time_re = re.compile(r"SimTimeSecs: ([0-9.]+)")
full_time_re = re.compile(r"FullTimeSecs: ([0-9.]+)")

# Begin output CSV
with open('performance_results.csv', 'w', newline='') as csv_file:
    csv_out = csv.writer(csv_file, delimiter=',')
    csv_out.writerow(("N_CELLS", "ECM_N", "SimTimeSecs", "FullTimeSecs"))

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
            sim_time_result = re.search(sim_time_re, output)
            full_time_result = re.search(full_time_re, output)
            if sim_time_result is None or full_time_result is None:
                print("N_CELLS %d, ECM_N %d output did not contain performance, skipping."%(N_CELLS, ECM_N))
                continue
            sim_time = float(sim_time_result.group(1))
            full_time = float(full_time_result.group(1))
            csv_out.writerow((N_CELLS, ECM_N, sim_time, full_time))

# Cleanup
if os.path.exists(temp_model_path):
    os.remove(temp_model_path)