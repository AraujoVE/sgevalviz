#!/bin/bash

bump_type="$1"

if [[ "$(uname)" == "Darwin" ]]; then
    sed_inplace=(-i '')
else
    sed_inplace=(-i)
fi

if [[ "$bump_type" == "major" || "$bump_type" == "minor" || "$bump_type" == "patch" || "$bump_type" == "rc" || "$bump_type" == "official" ]]; then
    file="pyproject.toml"
    file_content=$(grep -E '^version = ' "$file") || {
        echo "Could not find version line in $file"
        exit 1
    }
    official_regex="version = \"([0-9]+)\.([0-9]+)\.([0-9]+)\""
    test_regex="version = \"([0-9]+)\.([0-9]+)\.([0-9]+)rc([0-9]+)\""

    if [[ $file_content =~ $official_regex ]]; then
        major="${BASH_REMATCH[1]}"
        minor="${BASH_REMATCH[2]}"
        patch="${BASH_REMATCH[3]}"
 
        if [[ "$bump_type" == "rc" || "$bump_type" == "official" ]]; then
            echo "The file you are trying to update was previously an official version, you can bump it only to a new major, minor or patch version, and when you do this it will become a rc1 version of the new upgraded version. For example, upgrading a minor version of 3.2.11 will create 3.3.0rc1"
            exit 1
        else
            if [[ "$bump_type" == "major" ]]; then
                ((major++))
                minor=0
                patch=0
            elif [[ "$bump_type" == "minor" ]]; then
                ((minor++))
                patch=0
            else
                ((patch++))
            fi

            final_string="${major}.${minor}.${patch}rc1"
            sed "${sed_inplace[@]}" "s/version = \".*\"/version = \"${final_string}\"/" "$file"
        fi
    elif [[ $file_content =~ $test_regex ]]; then
        major="${BASH_REMATCH[1]}"
        minor="${BASH_REMATCH[2]}"
        patch="${BASH_REMATCH[3]}"
        rc="${BASH_REMATCH[4]}"
        
        if [[ "$bump_type" == "major" || "$bump_type" == "minor" || "$bump_type" == "patch" ]]; then
            echo "The file you are trying to update was previously an rc version, you can only bump the rc value or make it official to deploy. For example, upgrading a 3.2.11rc3 version with rc will make it 3.2.11rc4, while upgrading it with official will make it a 3.2.11 version"
            exit 1
        else               
            if [[ "$bump_type" == "rc" ]]; then
                ((rc++))
                final_string="${major}.${minor}.${patch}rc${rc}"
            else
                final_string="${major}.${minor}.${patch}"
            fi
            sed "${sed_inplace[@]}" "s/version = \".*\"/version = \"${final_string}\"/" "$file"
        fi
    else
        echo "file pyproject.toml doesn't have the correct version pattern"
        exit 1
    fi
else
    echo "The first argument must be either 'major', 'minor' or 'patch', if the current version is official. If the current version has 'rc', the first argument must be either 'rc' or 'official'."
    exit 1
fi

version_name="v${final_string}"
git add pyproject.toml
git commit -m "Bumping version to ${version_name}"
git push origin main
git tag -a "$version_name" -m "Release ${version_name}"
git push origin "$version_name"