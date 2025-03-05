#!/usr/bin/env python3

import json
import re
import subprocess

# Path to your devcontainer.json file
devcontainer_path = ".devcontainer/devcontainer.json"


def remove_comments(json_like_str):
    """Remove single line comments from JSON-like strings."""
    return re.sub(r"//.*$", "", json_like_str, flags=re.MULTILINE)


def check_code_command():
    """Check if 'code' or 'codium' command is available."""
    try:
        subprocess.run(
            ["code", "--version"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        return "code"
    except FileNotFoundError:
        try:
            subprocess.run(
                ["cursor", "--version"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return "cursor"
        except FileNotFoundError:
            print("Neither 'code' nor 'cursor' command is available.")
            return None


def install_vscode_extensions(devcontainer_path):
    # Check for available command
    code_cmd = check_code_command()
    if not code_cmd:
        return

    # Check if the devcontainer.json file exists
    try:
        with open(devcontainer_path, "r") as file:
            json_content = file.read()
            cleaned_json_content = remove_comments(json_content)
            devcontainer_config = json.loads(cleaned_json_content)
    except FileNotFoundError:
        print(f"The file {devcontainer_path} does not exist.")
        return
    except json.JSONDecodeError:
        print("Error decoding the JSON file.")
        return

    # Extract the list of extensions
    try:
        extensions = devcontainer_config["customizations"]["vscode"]["extensions"]
    except KeyError:
        print("No extensions found in the devcontainer.json file.")
        return

    if not extensions:
        print("No extensions found to install.")
        return

    # Install each extension using the VS Code/Cursor CLI
    for extension in extensions:
        try:
            # Skip GitHub extensions if using cursor
            # if code_cmd == "cursor" and extension.lower().startswith("github.copilot"):
            #     print(f"Skipping GitHub extension for cursor: {extension}")
            #     continue

            print(f"Installing extension: {extension}")
            subprocess.run([code_cmd, "--install-extension", extension], check=True)
        except subprocess.CalledProcessError:
            print(f"Failed to install extension: {extension}")


if __name__ == "__main__":
    install_vscode_extensions(devcontainer_path)
