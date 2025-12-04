import json
import subprocess
import sys
import os

PARAM_FILE = 'parameters 6 bolts.json'
MAIN_SCRIPT = 'main.py'

def read_params():
    with open(PARAM_FILE, 'r') as f:
        return json.load(f)

def write_params(params):
    with open(PARAM_FILE, 'w') as f:
        json.dump(params, f, indent=4)

def run_simulation():
    # Run main.py and capture output
    result = subprocess.run([sys.executable, MAIN_SCRIPT], capture_output=True, text=True)
    if result.returncode != 0:
        print("Error running main.py")
        print(result.stderr)
        return None, -999
    return parse_output(result.stdout)

def parse_output(output):
    lines = output.splitlines()
    ms_values = {
        't2_bear': [],
        't2_bear_th': [],
        't2_pull': [],
        't3_bear_th': [],
        't3_pull': []
    }
    
    start_parsing = False
    for line in lines:
        if "ID" in line and "X (m)" in line:
            start_parsing = True
            continue
        if start_parsing and line.strip().startswith("-"):
            continue
        if start_parsing:
            parts = line.split()
            if len(parts) >= 8:
                try:
                    # Columns: ID, X, Z, MS Bear(t2), MS Bear(t2 Th), MS Pull(t2), MS Bear(t3 Th), MS Pull(t3)
                    ms_values['t2_bear'].append(float(parts[3]))
                    ms_values['t2_bear_th'].append(float(parts[4]))
                    ms_values['t2_pull'].append(float(parts[5]))
                    ms_values['t3_bear_th'].append(float(parts[6]))
                    ms_values['t3_pull'].append(float(parts[7]))
                except ValueError:
                    pass
    
    # Flatten all MS values to check global constraint
    all_ms = []
    for v in ms_values.values():
        all_ms.extend(v)
    
    return ms_values, min(all_ms) if all_ms else -999

def optimize_parameter(param_keys, target_ms_keys, tolerance=1e-5):
    print(f"Optimizing {param_keys[-1]}...")
    
    # Initial bounds
    params = read_params()
    current_val = params
    for key in param_keys[:-1]:
        current_val = current_val[key]
    current_val = current_val[param_keys[-1]]
    
    low = 1e-5 # Minimum physical dimension
    high = current_val * 2 # Start with a reasonable upper bound
    
    # Check if we need to expand high
    # (Not implemented for now, assuming current is feasible or close)
    
    best_val = current_val
    
    for i in range(100): # Max iterations
        mid = (low + high) / 2
        
        # Update params
        params = read_params()
        curr = params
        for key in param_keys[:-1]:
            curr = curr[key]
        curr[param_keys[-1]] = mid
        write_params(params)
        
        ms_data, min_global_ms = run_simulation()
        
        if ms_data is None:
            print("Simulation failed.")
            break
            
        # Calculate min of TARGET MS
        target_vals = []
        for key in target_ms_keys:
            target_vals.extend(ms_data[key])
        
        min_target = min(target_vals) if target_vals else -999
        
        print(f"  Val: {mid:.6f}, Min Target MS: {min_target:.4f}, Global Min MS: {min_global_ms:.4f}")
        
        # Logic:
        # We want to minimize the parameter (mid).
        # Constraint: min_global_ms >= 0.
        # If min_global_ms < 0, we are too low. Increase low.
        # If min_global_ms >= 0, we are feasible. We can try to go lower.
        
        if min_target < 0:
            low = mid
        else:
            best_val = mid
            high = mid
            
            # If we are very close to 0 on the target, we might stop?
            # But we want to minimize the parameter as much as possible.
            # The limit is when min_global_ms = 0.
            # If min_target is the limiting factor, then min_target will be 0.
            # If something else is limiting, min_target might be > 0.
            # But the user said "change t3 until t3 values... are 0".
            # This implies t3 IS the limiting factor for t3 optimization.
            
            if abs(min_target) < tolerance:
                 break
                
    # Set to best found value
    params = read_params()
    curr = params
    for key in param_keys[:-1]:
        curr = curr[key]
    curr[param_keys[-1]] = best_val
    write_params(params)
    print(f"Set {param_keys[-1]} to {best_val}")

def main():
    # 1. Optimize t3
    # Target: t3_bear_th, t3_pull
    optimize_parameter(['geometry', 't3'], ['t3_bear_th', 't3_pull'])
    
    # 2. Optimize t2
    # Target: t2_bear, t2_bear_th, t2_pull
    optimize_parameter(['geometry', 't2'], ['t2_bear', 't2_bear_th', 't2_pull'])
    
    # 3. Optimize D_in
    # Target: All affected. Since D_in affects pull through and thermal, check all.
    # "change D_in until its closes to 0" -> minimize D_in until MS=0.
    optimize_parameter(['geometry', 'D_in'], ['t2_bear', 't2_bear_th', 't2_pull', 't3_bear_th', 't3_pull'])

if __name__ == "__main__":
    main()
