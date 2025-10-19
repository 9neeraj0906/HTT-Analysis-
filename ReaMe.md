# HTT Analysis CMS Open Source Data (2016)

This repository contains my CMS HTT analysis code in Python.

I have created a Docker container with 3 images:

* CMSSW_10_6_30
* Python
* ROOT

Data can be downloaded from the CMS Open Data portal based on the needs.

# Installation

Instructions and images are available at: [https://opendata.cern.ch/docs/cms-guide-docker](https://opendata.cern.ch/docs/cms-guide-docker)

Steps:

1. Download the images.
2. Create necessary project folders (Ideally separate directories for codes, datasets, results).
3. Follow the instructions from the above-mentioned page (Do not forget to mount the volumes i.e., the project directories).
4. Docker container is built.
5. Files inside the mounted volumes can now be accessed inside Docker.
6. NOTE: Disable Conda if you have it on your PC.

# Useful commands inside Docker

* `docker start -i <Docker name>` // Starts the Docker container
* `exit` // Exit Docker
* `cd /pythonn` // Navigates inside the 'pythonn' directory

# Git Integration

1. Go to the folder you want to integrate.
2. Initialize the repository:

   ```
   git init
   ```
3. Add files to the staging area:

   ```
   git add <filename>  # or git add . to add all files
   ```
4. Commit changes:

   ```
   git commit -m "Initial commit"
   ```
5. Add a remote repository:

   ```
   git remote add origin <repository-URL>
   ```
6. Push to a branch (e.g., main, temp, cppWay):

   ```
   git push -u origin <branch-name>
   ```

# Branch Management

* Create a new branch:

  ```
  git checkout -b <branch-name>
  ```
* Switch to an existing branch:

  ```
  git checkout <branch-name>
  ```
* Merge a branch into another (e.g., temp into main):

  ```
  git checkout main
  git merge temp
  ```

# Removing Files from a Branch

* To remove a file from a specific branch without affecting others:

  ```
  git checkout <branch-to-modify>
  git rm <filename>
  git commit -m "Remove file from this branch"
  git push origin <branch-to-modify>
  ```

# Git Authentication

* Use Personal Access Token (PAT) instead of password for HTTPS:

  ```
  git config --global credential.helper store  # Stores credentials for future use
  ```
* After configuring, you won't need to enter your PAT every time for this machine.

# Additional Tips

* `git status` shows which files are untracked, modified, or staged.
* `git log --oneline` gives a quick overview of commits.
* Always `git pull` before pushing to avoid conflicts:

  ```
  git pull origin <branch-name>
  ```
* Be careful when merging branches to avoid accidentally adding unwanted files.
* Personal workflow: Maintain `main` for finalized code, `temp` and `cppWay` for development branches.
* Only push files to branches where they should exist.
