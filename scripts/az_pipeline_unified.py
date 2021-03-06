import sys
import subprocess
import os
import argparse
import time

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run the AZ mapping '
                                     'pipeline with disambigutation')
    # read in the config file and perform initial setup
    parser.add_argument('mouse_mapping',
                        help="YAML file for mapping to mouse.")
    parser.add_argument('human_mapping',
                        help="YAML file for mapping to human.")
    parser.add_argument('disambiguate',
                        help="YAML file for disambiguation.")
    parser.add_argument('mouse_quantitation',
                        help="YAML file for running the mouse quantitation.")
    parser.add_argument('human_quantitation',
                        help="YAML file for running the human quantitation.")

    args = parser.parse_args()
    this_path = os.path.abspath(os.path.dirname(__file__))
    mapping_script_path = os.path.join(this_path, "mapping.py")
    mouse_mapping_path = os.path.abspath(args.mouse_mapping)
    human_mapping_path = os.path.abspath(args.human_mapping)
    mouse_mapping_cmd = ["python", mapping_script_path, mouse_mapping_path]
    human_mapping_cmd = ["python", mapping_script_path, human_mapping_path]
    procs = []
    for cmd in [mouse_mapping_cmd, human_mapping_cmd]:
        procs.append(subprocess.Popen(cmd))
    codes = [p.wait() for p in procs]
    if not all([c == 0 for c in codes]):
        print "One of the mapping scripts did not complete properly. Exiting."
        sys.exit(1)

    procs = []
    disambiguation_script_path = os.path.join(this_path, "disambiguate.py")
    disambiguation_path = os.path.abspath(args.disambiguate)
    disambiguation_cmd = ["python", disambiguation_script_path,
                          disambiguation_path]
    for cmd in [disambiguation_cmd]:
        procs.append(subprocess.Popen(cmd))
    codes = [p.wait() for p in procs]
    if not all([c == 0 for c in codes]):
        print "The disambiguation script did not complete properly. Exiting."
        sys.exit(1)

    procs = []
    quantitation_script_path = os.path.join(this_path, "quantitation.py")
    mouse_quantitation_path = os.path.abspath(args.mouse_quantitation)
    human_quantitation_path = os.path.abspath(args.human_quantitation)
    mouse_quantitation_cmd = ["python", quantitation_script_path, mouse_quantitation_path]
    human_quantitation_cmd = ["python", quantitation_script_path, human_quantitation_path]
    for cmd in [mouse_quantitation_cmd, human_quantitation_cmd]:
        procs.append(subprocess.Popen(cmd))
    codes = [p.wait() for p in procs]
    if not all([c == 0 for c in codes]):
        print "One of the quantitation scripts did not complete properly. Exiting."
        sys.exit(1)

    print "Run complete."
