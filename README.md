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

# Instruction to run the pipeline 
1. Run the server and give the qc method name which you want to run like hwe , mind , geno etc. (Any order of method is supported)
python Flower/qcServer.py --methods hwe snp_missingness ind_missingness maf_filter mind geno calculate_maf
2. Run the server and give the qc method name which you want to run like hwe , mind , geno etc. (Any order of method is supported)
python Flower/qcClient.py --client_id 1 --bed_file /Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bed --bim_file /Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.bim --fam_file /Users/sonamrathod/Documents/Rutgers/Project/Git_Repo/cc.qc3.fam --methods hwe snp_missingness ind_missingness maf_filter mind geno calculate_maf