from MKmodel.mk import load_microdosimetry_from_root, compute_basic_microdosimetry
from pathlib import Path
import os

# tree_name =Path("single_events")
# tree_name =Path("alpha_events")
tree_name =Path("events")

def get_all_file_path(folder_name):
    file_paths = []
    file_simple_path = []
    for root, dirs, files in os.walk(folder_name):
        dirs.sort()
        files.sort()
        for file in files:
            if file.endswith(".root"):
                file_paths.append(os.path.join(root, file))
                file_simple_path.append(file)
    return file_paths, file_simple_path

def draw_table(target_folder, tree_name):
    print("="*103)
    print(f"target_folder = {target_folder}, tree_name = {tree_name}")
    print("="*103)
    print(f"| {'target_file':<42} | {'z_d_F':<12} | {'z_d_D':<12} | {'z_n_F':<12} | {'z_n_D':<12} |")
    print("-"*107)
    file_paths, file_simple_paths = get_all_file_path(target_folder)
    for root_file, file_simple_path in zip(file_paths, file_simple_paths):
        if tree_name == Path("alpha_events") and "lu177" in root_file.lower():
            z_d, z_n, w, n_hits, n_entries = load_microdosimetry_from_root(root_file, Path("beta_events"))
        else:
            z_d, z_n, w, n_hits, n_entries = load_microdosimetry_from_root(root_file, tree_name)
        z_d_F, z_d_D, z_n_F, z_n_D = compute_basic_microdosimetry(z_d, z_n, w, n_hits)
        print(f"| {file_simple_path.split('.')[0]:<42} | {z_d_F:<12.6f} | {z_d_D:<12.6f} | {z_n_F:<12.6f} | {z_n_D:<12.6f} |")
        print("-"*107)

target_folders = [
    "./data/PC-ac225-phy-decay",
    "./data/V79-ac225-single-decay",
]

tree_names = [
    Path("single_events"),
    Path("alpha_events"),
    Path("events")
]

for target_folder in target_folders:
    for tree_name in tree_names:
        draw_table(target_folder, tree_name)
        print("\n\n")


    