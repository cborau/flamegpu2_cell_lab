from subprocess import check_output, CalledProcessError, STDOUT
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
init_time_re = re.compile(r"InitTimeSecs: ([0-9.]+)")
rtc_time_re = re.compile(r"RTCTimeSecs: ([0-9.]+)")
exit_time_re = re.compile(r"ExitTimeSecs: ([0-9.]+)")
step_avg_time_re = re.compile(r"StepTimeAvgSecs: ([0-9.]+)")

def run_model(N_CELLS, ECM_N):
    # Perform replacement of the two variables
    # @note These replacements are unstable if the template changes
    model_src = model_template.replace("<N_CELLS>", "%d"%(N_CELLS), 1).replace("<ECM_N>", "%d"%(ECM_N), 1)
    with open(temp_model_path, "w") as model_file:
        model_file.write(model_src)
    # Execute (might require python3)
    try:
        return check_output(["python", temp_model_path], stderr=STDOUT).decode()
    except CalledProcessError as e:
        return e.output.decode()

# Perform an empty run first, to force RTC to compile
run_model(100, 10)

# Begin output CSV
with open('performance_results.csv', 'w', newline='') as csv_file:
    csv_out = csv.writer(csv_file, delimiter=',')
    csv_out.writerow(("N_CELLS", "ECM_N", "SimTimeSecs", "FullTimeSecs", "InitTimeSecs", "RTCInitTimeSecs", "ExitTimeSecs", "StepTimeAvgSecs"))

    # For each configuration
    N_CELLS_sweep = [1000, 10000, 100000, 1000000]
    ECM_N_sweep = range(10,101,10)
    for N_CELLS in N_CELLS_sweep:
        for ECM_N in ECM_N_sweep:
            print("N_CELLS %d, ECM_N %d"%(N_CELLS, ECM_N), end='')
            # Perform replacement of the two variables
            # @note These replacements are unstable if the template changes
            output = run_model(N_CELLS, ECM_N)
            # Collect/validate performance results
            sim_time_result = re.search(sim_time_re, output)
            full_time_result = re.search(full_time_re, output)
            init_time_result = re.search(init_time_re, output)
            rtc_time_result = re.search(rtc_time_re, output)
            exit_time_result = re.search(exit_time_re, output)
            step_avg_time_result = re.search(step_avg_time_re, output)
            if (sim_time_result is None or 
                full_time_result is None or 
                init_time_result is None or 
                rtc_time_result is None or 
                exit_time_result is None or 
                step_avg_time_result is None):
                print(", output did not contain performance, skipping.")
                continue
            sim_time = float(sim_time_result.group(1))
            full_time = float(full_time_result.group(1))
            init_time = float(init_time_result.group(1))
            rtc_time = float(rtc_time_result.group(1))
            exit_time = float(exit_time_result.group(1))
            step_avg_time = float(step_avg_time_result.group(1))
            csv_out.writerow((N_CELLS, ECM_N, sim_time, full_time, init_time, rtc_time, exit_time, step_avg_time))
            csv_file.flush()
            print(", sim %fs, full %fs, rtc %fs, init %fs, step avg %s"%(sim_time, full_time, rtc_time, init_time, step_avg_time))

# Cleanup
if os.path.exists(temp_model_path):
    os.remove(temp_model_path)