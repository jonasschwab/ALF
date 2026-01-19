#!/usr/bin/env python3
import argparse
import json
import os
import yaml

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare test directories for test vs ED based on specifications.")
    parser.add_argument('--specs_file', default="testsuite/test_vs_ed/test_specs.yaml", help="Path to the test specifications YAML file.")
    parser.add_argument('--output_dir', default="prepared-test-directories", help="Directory to store the prepared test directories.")
    args = parser.parse_args()

    with open(args.specs_file, 'r', encoding='UTF-8') as f:
        test_specs = yaml.safe_load(f)
    
    simulation_matrix = []

    for test_name, test_spec in test_specs.items():
        for env_name, env_spec in test_spec['environments'].items():
            test_dir = f"{args.output_dir}/{test_name}_{env_name}"
            os.makedirs(test_dir)
            with open(os.path.join(test_dir, "spec.yaml"), 'w', encoding='UTF-8') as f:
                f.write(yaml.dump({
                    'test_name': test_name,
                    'env_name': env_name,
                    **test_spec,
                    }))
            for i, _ in enumerate(test_spec["sim_dicts"]):
                os.makedirs(os.path.join(test_dir, f"{i+1}"))
                with open(os.path.join(test_dir, f"{i+1}", "spec.yaml"), 'w', encoding='UTF-8') as f:
                    f.write(yaml.dump({
                        'machine': env_spec['variables']['MACHINE'],
                        'ham_name': test_spec["ham_name"],
                        'sim_dict': test_spec["sim_dicts"][i],
                    }))
                simulation_matrix.append({
                    'test_name': test_name,
                    'env_name': env_name,
                    'machine': env_spec['variables']['MACHINE'],
                    'image': env_spec['image'],
                    'CI_NODE_INDEX': i+1,
                })
    
    with open(os.getenv('GITHUB_OUTPUT', 'test.json'), 'w+', encoding='UTF-8') as f:
      f.write(f'simulation_matrix={json.dumps(simulation_matrix)}\n')
