#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import re

# Try importing yaml, fallback to basic parsing if not available (though PyYAML should be present)
try:
    import yaml
except ImportError:
    yaml = None

# Mappings of config keys to package details
PACKAGE_MAPPING = {
    "datamodel": {
        "package_name": "pfs-datamodel",
        "repo_url": "https://github.com/Subaru-PFS/datamodel.git",
        "required": True
    },
    "pfs_utils": {
        "package_name": "pfs-utils",
        "repo_url": "https://github.com/Subaru-PFS/pfs_utils.git",
        "required": True
    },
    "ics_cobraCharmer": {
        "package_name": "ics-cobraCharmer",
        "repo_url": "https://github.com/Subaru-PFS/ics_cobraCharmer.git",
        "required": False
    },
    "ics_cobraOps": {
        "package_name": "ics-cobraOps",
        "repo_url": "https://github.com/Subaru-PFS/ics_cobraOps.git",
        "required": True
    },
    "ets_fiberalloc": {
        "package_name": "ets-fiber-assigner",
        "repo_url": "https://github.com/Subaru-PFS/ets_fiberalloc.git",
        "required": True
    },
    "pfs_instdata": {
        "package_name": "pfs-instdata",
        "repo_url": "https://github.com/Subaru-PFS/pfs_instdata.git",
        "required": True
    },
    "pfs_design_tool": {
        "package_name": "pfs_design_tool",
        "repo_url": "https://github.com/Subaru-PFS/ets_pointing.git",
        "required": True
    },
    "pfs_obsproc_planning_tools": {
        "package_name": "pfs_obsproc_planning",
        "repo_url": "https://github.com/Subaru-PFS/pfs_obsproc_planning_tools.git",
        "required": True,
        "no_deps": True  # qplan depends on ginga>=6.1.0 which conflicts with other packages
    }
}

def apply_validation_patch(venv_dir):
    """Applies a bugfix patch to validation.py to correctly handle missing fluxes."""
    import glob
    search_pattern = os.path.join(venv_dir, "lib", "python3.*", "site-packages", "pfs_obsproc_planning", "utils", "validation.py")
    matches = glob.glob(search_pattern)
    if not matches:
        print("Warning: Could not find validation.py to apply patch.")
        return
        
    validation_file = matches[0]
    print(f"\nApplying validation.py bugfix patch to: {validation_file}")
    
    try:
        with open(validation_file, "r") as f:
            content = f.read()
            
        old_block = '''        elif t == 1:  # TargetType.SCIENCE
            if filt is not None and filt in fl:
                # Preserve the original branching behaviour (note: original used `if not np.nan`)
                pfsflux_l.append(
                    a[fl.index(filt)]
                    if not np.isnan(a[fl.index(filt)])
                    else b[fl.index(filt)]
                )
                pfsflux_f.append(filt)
            else:
                indices = [i for i, item in enumerate(fl) if ~np.isin(item, ["none", "nan"])]
                pfsflux_l.append(
                    a[indices[0]] if not np.isnan(a[indices[0]]) else b[indices[0]]
                )
                pfsflux_f.append(fl[indices[0]])'''

        new_block = '''        elif t == 1:  # TargetType.SCIENCE
            filt_idx = fl.index(filt) if (filt is not None and filt in fl) else -1
            if filt_idx >= 0 and not (np.isnan(a[filt_idx]) and np.isnan(b[filt_idx])):
                pfsflux_l.append(
                    a[filt_idx] if not np.isnan(a[filt_idx]) else b[filt_idx]
                )
                pfsflux_f.append(filt)
            else:
                indices = [
                    i for i, item in enumerate(fl)
                    if not np.isin(item, ["none", "nan"]) and not (np.isnan(a[i]) and np.isnan(b[i]))
                ]
                if len(indices) > 0:
                    idx = indices[0]
                    pfsflux_l.append(
                        a[idx] if not np.isnan(a[idx]) else b[idx]
                    )
                    pfsflux_f.append(fl[idx])
                else:
                    pfsflux_l.append(np.nan)
                    pfsflux_f.append("none")'''
                    
        if old_block in content:
            new_content = content.replace(old_block, new_block)
            with open(validation_file, "w") as f:
                f.write(new_content)
            print("Successfully applied patch to validation.py")
        else:
            print("Warning: Could not apply patch, target code block not found in validation.py. It might have been already patched or updated upstream.")
    except Exception as e:
        print(f"Error applying patch to validation.py: {e}")

def apply_plot_pfsdesign_patch(venv_dir):
    """Applies a patch to plot_pfsDesign.py to remove the newline inside the Std star label."""
    import glob
    search_pattern = os.path.join(venv_dir, "lib", "python3.*", "site-packages", "pfs_obsproc_planning", "utils", "plot_pfsDesign.py")
    matches = glob.glob(search_pattern)
    if not matches:
        print("Warning: Could not find plot_pfsDesign.py to apply patch.")
        return
        
    plot_file = matches[0]
    print(f"\nApplying plot_pfsDesign.py patch to: {plot_file}")
    
    try:
        with open(plot_file, "r") as f:
            content = f.read()
            
        old_block = '        label=f"Std star ({filtername_std}, {len(std)}, \\n {proposal_detail})",'
        new_block = '        label=f"Std star ({filtername_std}, {len(std)}, {proposal_detail})",'
        
        if old_block in content:
            new_content = content.replace(old_block, new_block)
            with open(plot_file, "w") as f:
                f.write(new_content)
            print("Successfully applied patch to plot_pfsDesign.py")
        else:
            print("Warning: Could not apply patch to plot_pfsDesign.py, target label block not found.")
    except Exception as e:
        print(f"Error applying patch to plot_pfsDesign.py: {e}")

def parse_yaml_fallback(file_path):
    """Fallback parser if PyYAML is not installed."""
    pfs_deps = {}
    gurobi_ver = None
    
    in_pfs = False
    in_gurobi = False
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Check sections
            if line.startswith('pfs:'):
                in_pfs = True
                in_gurobi = False
                continue
            elif line.startswith('gurobi:'):
                in_gurobi = True
                in_pfs = False
                continue
            elif line.endswith(':') or (':' in line and line.split(':')[0].strip() not in pfs_deps and not in_pfs and not in_gurobi):
                # Another top-level section
                in_pfs = False
                in_gurobi = False
            
            # Parse keys
            if in_pfs:
                # Expecting indented key-value pairs
                m = re.match(r"^([a-zA-Z0-9_]+)\s*:\s*(.+)$", line)
                if m:
                    key = m.group(1).strip()
                    val = m.group(2).split('#')[0].strip().strip('"').strip("'")
                    if val.lower() == 'null' or val.lower() == 'none':
                        val = None
                    pfs_deps[key] = val
            elif in_gurobi:
                m = re.match(r"^version\s*:\s*(.+)$", line)
                if m:
                    gurobi_ver = m.group(1).split('#')[0].strip().strip('"').strip("'")
                    
    return pfs_deps, gurobi_ver

def load_config(file_path):
    """Loads configuration from YAML file."""
    if not os.path.exists(file_path):
        print(f"Error: Configuration file not found at {file_path}", file=sys.stderr)
        sys.exit(1)
        
    pfs_deps = {}
    gurobi_ver = None
    
    if yaml is not None:
        try:
            with open(file_path, 'r') as f:
                config = yaml.safe_load(f)
            if 'pfs' in config:
                pfs_deps = config['pfs']
            if 'gurobi' in config and 'version' in config['gurobi']:
                gurobi_ver = config['gurobi']['version']
        except Exception as e:
            print(f"Warning: PyYAML failed to parse file ({e}). Falling back to manual parser.", file=sys.stderr)
            pfs_deps, gurobi_ver = parse_yaml_fallback(file_path)
    else:
        pfs_deps, gurobi_ver = parse_yaml_fallback(file_path)
        
    return pfs_deps, gurobi_ver

def parse_week_tag(tag):
    """Parses a weekly tag (e.g. w.2026.22) into (year, week) tuple."""
    if tag.startswith("w."):
        parts = tag[2:].split(".")
        if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
            return int(parts[0]), int(parts[1])
    return None

def parse_semver_tag(tag):
    """Parses a semantic version tag (e.g. v1.8.72) into numerical tuple."""
    t = tag.lstrip("vV")
    parts = t.split(".")
    numeric_parts = []
    for p in parts:
        if p.isdigit():
            numeric_parts.append(int(p))
        else:
            m = re.match(r"^(\d+)", p)
            if m:
                numeric_parts.append(int(m.group(1)))
            else:
                break
    if numeric_parts:
        return tuple(numeric_parts)
    return None

def tag_sort_key(tag):
    """Sorting key for tags, handling weekly and semver tags appropriately."""
    wk = parse_week_tag(tag)
    if wk is not None:
        return (0, wk) # Weekly tags first
    sv = parse_semver_tag(tag)
    if sv is not None:
        return (1, sv) # Semver tags second
    return (-1, tag) # Fallback lexicographical

def parse_constraint(constraint_str):
    """Parses operator and version from a constraint string."""
    constraint_str = constraint_str.strip()
    for op in [">=", "<=", ">", "<", "=="]:
        if constraint_str.startswith(op):
            return op, constraint_str[len(op):].strip()
    # If no operator is found, assume direct equality (exact tag/branch/commit)
    return "==", constraint_str

def match_version(tag, op, target_ver_str):
    """Checks if a tag satisfies a given operator and target version."""
    is_weekly = target_ver_str.startswith("w.")
    if is_weekly:
        tag_parsed = parse_week_tag(tag)
        target_parsed = parse_week_tag(target_ver_str)
        if tag_parsed is None or target_parsed is None:
            return False
    else:
        tag_parsed = parse_semver_tag(tag)
        target_parsed = parse_semver_tag(target_ver_str)
        if tag_parsed is None or target_parsed is None:
            return False
            
    if op == "==":
        return tag_parsed == target_parsed
    elif op == ">":
        return tag_parsed > target_parsed
    elif op == ">=":
        return tag_parsed >= target_parsed
    elif op == "<":
        return tag_parsed < target_parsed
    elif op == "<=":
        return tag_parsed <= target_parsed
    return False

def get_git_tags(repo_url):
    """Queries the remote Git repository for all tags."""
    cmd = ["git", "ls-remote", "--tags", repo_url]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error querying remote tags for {repo_url}: {e.stderr}", file=sys.stderr)
        return []
        
    tags = []
    for line in res.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        ref = parts[1]
        if ref.startswith("refs/tags/"):
            tag = ref[len("refs/tags/"):]
            if tag.endswith("^{}"):
                tag = tag[:-3]
            if tag not in tags:
                tags.append(tag)
    return tags

def resolve_git_ref(package_key, spec, repo_url):
    """Resolves the version specification to a specific Git tag or branch."""
    if not spec:
        return None
        
    op, target_ver = parse_constraint(spec)
    
    # If the operator is '==' and target_ver doesn't look like a constraint that needs resolution
    # (or if it's a specific tag/branch), we can check remote tags but we can also just use it directly.
    if op == "==":
        # Check if the tag starts with 'v' or is a weekly tag
        return target_ver
        
    print(f"Resolving constraint '{spec}' for {package_key}...")
    tags = get_git_tags(repo_url)
    if not tags:
        print(f"Warning: Could not fetch tags for {package_key}. Falling back to installing raw spec @{target_ver}.", file=sys.stderr)
        return target_ver
        
    matching_tags = [t for t in tags if match_version(t, op, target_ver)]
    if not matching_tags:
        print(f"Warning: No tags matched constraint '{spec}' for {package_key}. Falling back to default remote HEAD.", file=sys.stderr)
        return None
        
    # Sort tags to find the latest
    matching_tags.sort(key=tag_sort_key, reverse=True)
    resolved_tag = matching_tags[0]
    print(f"-> Resolved '{spec}' to tag '{resolved_tag}'")
    return resolved_tag

def cleanup_bin_directory(venv_dir):
    """Move non-essential executables in .venv/bin to an extras directory."""
    import shutil
    bin_dir = os.path.join(venv_dir, "bin")
    if not os.path.isdir(bin_dir):
        return
        
    extras_dir = os.path.join(bin_dir, "extras")
    os.makedirs(extras_dir, exist_ok=True)
    
    keep_list = [
        "python", "python3", "pip", "pip3", "uv", "activate", "activate.csh", "activate.fish", "Activate.ps1",
        "run_netflow", "make_pfs_design", "optimize_hex_fov", 
        "optimize_hex_fov_with_guidestars", "optimize_hex_fov_local_search",
        "merge_target_csv", "check_sky_sectors"
    ]
    
    moved_count = 0
    for filename in os.listdir(bin_dir):
        if filename == "extras":
            continue
            
        file_path = os.path.join(bin_dir, filename)
        
        # Keep if exactly matches the name, or if it matches python3.X, pip3.X, etc.
        should_keep = False
        if filename in keep_list:
            should_keep = True
        elif filename.startswith("python3.") or filename.startswith("pip3."):
            should_keep = True
            
        if not should_keep and os.path.isfile(file_path):
            shutil.move(file_path, os.path.join(extras_dir, filename))
            moved_count += 1
            
    if moved_count > 0:
        print(f"Cleaned up bin directory: moved {moved_count} extra executables to {extras_dir}")

def main():
    parser = argparse.ArgumentParser(description="Create a python virtual environment with uv based on PFS config dependencies.")
    parser.add_argument("config_file", help="Path to config.yaml (e.g. spt_ssp_observation/runs/2026-07/config.yaml)")
    parser.add_argument("-o", "--venv-dir", default=".venv", help="Path/directory to create the virtual environment (default: .venv)")
    parser.add_argument("-p", "--python", default="3.12", help="Python version to use for virtual environment (default: 3.12)")
    parser.add_argument("--use-local", action="store_true", help="Install local ets_fiberalloc and ics_cobraOps subdirectories in editable mode if present")
    parser.add_argument("--dry-run", action="store_true", help="Print actions and commands without executing them")
    
    args = parser.parse_args()
    
    pfs_deps, gurobi_ver = load_config(args.config_file)
    
    print("=" * 60)
    print(f"PFS Dependency Resolver for config: {args.config_file}")
    print("=" * 60)
    
    # 1. Create Virtual Environment Command
    venv_cmd = ["uv", "venv", args.venv_dir, "--python", args.python]
    print(f"Virtual environment command: {' '.join(venv_cmd)}")
    
    if not args.dry_run:
        print("Creating virtual environment...")
        subprocess.run(venv_cmd, check=True)
    
    # 2. Build install list
    install_targets = []
    
    # Check each package in package mapping
    for key, info in PACKAGE_MAPPING.items():
        spec = pfs_deps.get(key)
        if spec is None and not info.get("required", False):
            # Check case differences or ignore if null, and not required
            continue
            
        pkg_name = info["package_name"]
        repo_url = info["repo_url"]
        no_deps = info.get("no_deps", False)
        
        # Check if local installation is requested and available
        is_local_target = args.use_local and key in ["ets_fiberalloc", "ics_cobraOps"]
        if is_local_target and os.path.isdir(key):
            print(f"Using local directory for {key} (editable mode)")
            install_targets.append(("-e {key}", no_deps))
        else:
            resolved_ref = resolve_git_ref(key, spec, repo_url) if spec else None
            if resolved_ref:
                install_targets.append((f"{pkg_name} @ git+{repo_url}@{resolved_ref}", no_deps))
            else:
                install_targets.append((f"{pkg_name} @ git+{repo_url}", no_deps))
                
    # Add packages required for targetdb_calib.py
    print("Adding packages required for targetdb_calib.py...")
    targetdb_packages = [
        "sshtunnel",
        "paramiko<4.0.0",
        "sqlalchemy",
        "toml",
        "astropy",
        "pandas",
        "tqdm",
        "psycopg2-binary"
    ]
    for pkg in targetdb_packages:
        install_targets.append((pkg, False))

    # Add Gurobi if specified
    if gurobi_ver:
        # Check Python and Gurobi compatibility (gurobipy < 11 does not support Python >= 3.12)
        try:
            py_parts = [int(p) for p in args.python.split(".") if p.isdigit()]
            is_py312_or_newer = len(py_parts) >= 2 and (py_parts[0] > 3 or (py_parts[0] == 3 and py_parts[1] >= 12))
        except Exception:
            is_py312_or_newer = False
            
        gurobi_parts = parse_semver_tag(gurobi_ver)
        is_gurobi_old = gurobi_parts and len(gurobi_parts) >= 1 and gurobi_parts[0] < 11
        
        if is_py312_or_newer and is_gurobi_old:
            print(f"Warning: gurobipy=={gurobi_ver} does not support Python {args.python} (requires Python <= 3.11).")
            print(f"-> Automatically upgrading gurobipy requirement to >= 11.0.0 to support Python {args.python}.")
            install_targets.append(("gurobipy>=11.0.0", False))
        else:
            print(f"Adding Gurobi version: {gurobi_ver}")
            install_targets.append((f"gurobipy=={gurobi_ver}", False))
            
    # Install the current directory as an editable package so that pyproject.toml scripts are added to bin
    install_targets.append(("-e .", False))
        
    if not install_targets:
        print("No packages to install.")
        return

    # Split into normal and no-deps targets
    normal_targets = [t for t, nd in install_targets if not nd]
    nodeps_targets = [t for t, nd in install_targets if nd]
        
    # Generate overrides file to resolve potential nested VCS URL conflicts
    overrides_file = None
    git_targets = [t for t in normal_targets if " @ git+" in t]
    
    if git_targets:
        overrides_dir = os.path.abspath(args.venv_dir)
        if not args.dry_run:
            os.makedirs(overrides_dir, exist_ok=True)
        overrides_file = os.path.join(overrides_dir, "uv-overrides.txt")
        print(f"\nWriting overrides to: {overrides_file}")
        if not args.dry_run:
            with open(overrides_file, "w") as f:
                for target in git_targets:
                    f.write(f"{target}\n")
                    
    # Execute uv pip install
    # Note: we need to use the python/uv path in the newly created venv
    uv_executable = "uv"
    pip_cmd = [uv_executable, "pip", "install"]
    
    # Set the VIRTUAL_ENV environment variable to make sure uv uses the correct virtual env
    env = os.environ.copy()
    env["VIRTUAL_ENV"] = os.path.abspath(args.venv_dir)
    
    print("\nPackages to be installed:")
    for target in normal_targets:
        print(f"  - {target}")
    if nodeps_targets:
        print("Packages to be installed (--no-deps):")
        for target in nodeps_targets:
            print(f"  - {target}")
        
    # For sub-arguments containing space (like '-e path'), we should split or append correctly
    flat_targets = []
    for t in normal_targets:
        if t.startswith("-e "):
            flat_targets.extend(["-e", t[3:]])
        else:
            flat_targets.append(t)
            
    full_pip_cmd = pip_cmd + flat_targets
    if overrides_file:
        full_pip_cmd.extend(["--overrides", overrides_file])
        
    print(f"\nInstallation command:\n{' '.join(full_pip_cmd)}")
    
    if not args.dry_run:
        print("\nInstalling packages...")
        subprocess.run(full_pip_cmd, env=env, check=True)

    # Install no-deps packages separately
    if nodeps_targets:
        flat_nodeps = []
        for t in nodeps_targets:
            if t.startswith("-e "):
                flat_nodeps.extend(["-e", t[3:]])
            else:
                flat_nodeps.append(t)
        nodeps_pip_cmd = pip_cmd + ["--no-deps"] + flat_nodeps
        print(f"\nNo-deps installation command:\n{' '.join(nodeps_pip_cmd)}")
        if not args.dry_run:
            print("Installing no-deps packages...")
            subprocess.run(nodeps_pip_cmd, env=env, check=True)

    if not args.dry_run:
        apply_validation_patch(os.path.abspath(args.venv_dir))
        apply_plot_pfsdesign_patch(os.path.abspath(args.venv_dir))
        cleanup_bin_directory(os.path.abspath(args.venv_dir))
        print("\nVirtual environment setup completed successfully!")
    else:
        print("\nDry-run mode. No changes were made.")

if __name__ == "__main__":
    main()
