"""
git_setup.py — A simple helper to initialize git and commit all files in this project.
Run this script using Python if you want to automate staging and committing:
    python git_setup.py
"""

import subprocess
import sys
from pathlib import Path


def run_command(cmd, shell=True):
    print(f"Executing: {cmd}")
    res = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error executing command: {res.stderr.strip()}", file=sys.stderr)
        return False
    print(res.stdout.strip())
    return True


def main():
    proj_dir = Path(__file__).parent.resolve()
    print(f"Initializing Git repository in: {proj_dir}")

    # 1. git init
    if not run_command("git init"):
        sys.exit(1)

    # 2. git add .
    if not run_command("git add ."):
        sys.exit(1)

    # 3. git commit
    if not run_command('git commit -m "feat: initial sample QC implementation"'):
        print("\nNote: Commit failed. Please make sure your git user.name and user.email are configured:")
        print("  git config --global user.name \"Your Name\"")
        print("  git config --global user.email \"you@example.com\"")
        sys.exit(1)

    print("\n🎉 Git initialized and all files successfully committed!")
    print("\nNext steps to push to your GitHub account:")
    print("1. Create a new repository on GitHub named 'sample_QC'")
    print("2. Run the following commands in your terminal:")
    print("   git remote add origin https://github.com/<YOUR_GITHUB_USERNAME>/sample_QC.git")
    print("   git branch -M main")
    print("   git push -u origin main")


if __name__ == "__main__":
    main()
