import os
import subprocess

# Set path to Stride executable file
stride_path = './stride'

# Define function to run Stride and extract predicted dihedral angles and secondary structure
def get_stride_results(pdb_file):
    # Run Stride and capture output
    stride_output = subprocess.check_output([stride_path, pdb_file])
    # Convert output to string and split into lines
    stride_lines = stride_output.decode('utf-8').split('\n')
    # Extract predicted dihedral angles and secondary structure from output
    phi_psi_angles = {}
    secondary_structure = {}
    for line in stride_lines:
        if line.startswith('ASG'):
            res_id = line[7:11].strip()
            phi = line[21:29].strip()
            psi = line[33:41].strip()
            phi_psi_angles[res_id] = (phi, psi)
        elif line.startswith('STR'):
            res_id = line[7:11].strip()
            ss = line[24]
            secondary_structure[res_id] = ss
    return phi_psi_angles, secondary_structure

# Define function to calculate accuracy of predicted dihedral angles and secondary structure
def calculate_accuracy(true_angles, predicted_angles, true_ss, predicted_ss):
    # Calculate accuracy of dihedral angles
    phi_psi_accuracy = {}
    for res_id in true_angles:
        true_phi, true_psi = true_angles[res_id]
        predicted_phi, predicted_psi = predicted_angles.get(res_id, (None, None))
        if predicted_phi is None or predicted_psi is None:
            accuracy = 0.0
        else:
            phi_error = abs(float(true_phi) - float(predicted_phi))
            psi_error = abs(float(true_psi) - float(predicted_psi))
            accuracy = 1.0 - (phi_error + psi_error) / 360.0
        phi_psi_accuracy[res_id] = accuracy
    phi_psi_mean_accuracy = sum(phi_psi_accuracy.values()) / len(phi_psi_accuracy)
    # Calculate accuracy of secondary structure
    ss_accuracy = {}
    for res_id in true_ss:
        true_ss_code = true_ss[res_id]
        predicted_ss_code = predicted_ss.get(res_id, None)
        if predicted_ss_code is None:
            accuracy = 0.0
        else:
            accuracy = 1.0 if true_ss_code == predicted_ss_code else 0.0
        ss_accuracy[res_id] = accuracy
    ss_mean_accuracy = sum(ss_accuracy.values()) / len(ss_accuracy)
    return phi_psi_mean_accuracy, ss_mean_accuracy

pdb_file = 'protein.pdb'
predicted_phi_psi_angles, predicted_secondary_structure = get_stride_results(pdb_file)
phi_psi_accuracy, ss_accuracy = calculate_accuracy()
print(f"Phi/Psi accuracy: {phi_psi_accuracy:.3f}")
print(f"Secondary structure accuracy: {ss_accuracy:.3f}")
