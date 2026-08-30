import yaml
import sys
import os
from pathlib import Path

def key_aware_merge(dict1, dict2):
    for key, value in dict2.items():
        if isinstance(value, dict) and key in dict1 and isinstance(dict1[key], dict):
            key_aware_merge(dict1[key], value)
            
        elif isinstance(value, list) and key in dict1 and isinstance(dict1[key], list):
            # Map existing list dictionary top keys to their index in the list
            existing_indices = {}
            for idx, item in enumerate(dict1[key]):
                if isinstance(item, dict):
                    for sub_k in item.keys():
                        existing_indices[sub_k] = idx

            # Process new list items
            for item in value:
                if isinstance(item, dict):
                    matched = False
                    for sub_k, sub_v in item.items():
                        if sub_k in existing_indices:
                            # Key already exists: replace/merge the item in place
                            target_idx = existing_indices[sub_k]
                            if isinstance(dict1[key][target_idx][sub_k], dict) and isinstance(sub_v, dict):
                                key_aware_merge(dict1[key][target_idx][sub_k], sub_v)
                            else:
                                dict1[key][target_idx][sub_k] = sub_v
                            matched = True
                    if not matched:
                        dict1[key].append(item)
                else:
                    if item not in dict1[key]:
                        dict1[key].append(item)
        else:
            dict1[key] = value
    return dict1

if __name__ == "__main__":
    assert len(sys.argv)==2, f"Usage: {sys.argv[0]} gamsconfig_newsolver.yaml [> gamsconfig.yaml]"
    assert os.path.exists(sys.argv[1]), f"File {sys.argv[1]} does not exist!"
    if not os.path.exists("gamsconfig.yaml"):
        print(Path(sys.argv[1]).read_text(), end="")
    else:
        # Run the merge
        with open('gamsconfig.yaml', 'r') as f1, open(sys.argv[1], 'r') as f2:
            data1 = yaml.safe_load(f1) or {}
            data2 = yaml.safe_load(f2) or {}
        
        merged = key_aware_merge(data1, data2)
        yaml.dump(merged, sys.stdout, default_flow_style=False, sort_keys=False)