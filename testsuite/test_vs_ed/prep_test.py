#!/usr/bin/env python3
import json
import os
import sys
import yaml
import copy

if __name__ == "__main__":
    try:
        specs_file = sys.argv[1]
    except IndexError:
        specs_file = "testsuite/test_vs_ed/test_specs.yaml"
    with open(specs_file, 'r', encoding='UTF-8') as f:
        test_specs = yaml.safe_load(f)
    
    compile_matrix = []
    simulation_matrix = []
    analysis_matrix = []

    for test_name, test_spec in test_specs.items():
        for env_name, env_spec in test_spec['environments'].items():
            test_dir = f"testsuite/test_vs_ed/{test_name}_{env_name}"
            os.makedirs(test_dir)
            with open(os.path.join(test_dir, "spec.yaml"), 'w', encoding='UTF-8') as f:
                f.write(yaml.dump({
                    'test_name': test_name,
                    'env_name': env_name,
                    **test_spec,
                    }))
            compile_matrix.append({
                'test_name': test_name,
                'env_name': env_name,
                'machine': env_spec['variables']['MACHINE'],
                'image': env_spec['image'],
            })
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
      f.write(f'compile_matrix={json.dumps(compile_matrix)}\n')
      f.write(f'simulation_matrix={json.dumps(simulation_matrix)}')

    with open('generated-config.yml', 'w', encoding='UTF-8') as f:
        f.write(yaml.dump(pipeline_config))
