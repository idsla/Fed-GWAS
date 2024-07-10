Federated GWAS Pipeline

# Poetry Project setup

Install pipx in Windows
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
Invoke-RestMethod -Uri https://get.scoop.sh | Invoke-Expression

# Install pipx - C:\Users\xxx\.local\bin (path)
scoop install pipx
```

Install Poetry
```bash
# Install peotry - in git bash
pipx install poetry
pipx ensurepath
```

Install Python > 3.11.0 and follow the below workflow for using poetry

```bash
poetry completions bash >> ~/.bash_completion  # optional
poetry config virtualenvs.in-project true
poetry self update
# first time use - create virtual environment using poetry
poetry env use python

# activate enviroment
poetry shell

# install packages from poetry.lock file
poetry install

# add new packages
poetry add <package-name>

# add new pakcages to dev group - this group is for those packages using for helping our development, will not be included in production version
peory add <package_name> --group dev

# remove packages
poetry remove <package-name>

# deactivate environment
exit
```
