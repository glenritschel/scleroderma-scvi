# Full WSL + Paper1 Environment Bootstrap (Clean)

This document recreates the full development environment for the
`scleroderma-scvi` project on a fresh WSL Ubuntu install.

---

## System prerequisites (one time)

```bash
sudo apt update
sudo apt install -y \
    curl \
    ca-certificates \
    tar \
    bzip2 \
    build-essential \
    g++ \
    make
````

---

## Install micromamba

```bash
cd ~
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
sudo mv bin/micromamba /usr/local/bin/
micromamba shell init -s bash
source ~/.bashrc
```

Verify:

```bash
micromamba --version
```

---

## Fix broken conda config (if present)

```bash
mv ~/.condarc ~/.condarc.bak
```

---

## Generate SSH key (WSL)

```bash
ssh-keygen -t ed25519 -C "glen@DESKTOP-AM7VVAT"
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519
cat ~/.ssh/id_ed25519.pub
```

Add the public key to:

GitHub → Settings → SSH and GPG keys → New SSH key

Test:

```bash
ssh -T git@github.com
```

---

## Clone repos

```bash
mkdir -p ~/repos
cd ~/repos

git clone git@github.com:glenritschel/scleroderma-scvi.git
git clone git@github.com:glenritschel/ebv-immune-scarring.git
```

---

## Create Paper 1 environment

```bash
cd ~/repos/scleroderma-scvi
micromamba env create -n ssc-scvi -f environment.yml
micromamba activate ssc-scvi
```

---

## Register Jupyter kernel (recommended)

```bash
python -m pip install ipykernel
python -m ipykernel install --user --name ssc-scvi --display-name "ssc-scvi"
```

---

## Validate stack

```bash
python -c "import scanpy, scvi, anndata, annoy; print('environment OK')"
```

---

## Operational rules (important)

Always shut down WSL before reboot/sleep:

```bash
wsl --shutdown
```

Never unplug the external SSD while Windows is running.
Do not back up live `X:\wsl` VHDX files.
Keep active repos inside WSL (`~/repos`).

---

## Result

* Clean WSL on external SSD
* SSH GitHub access
* Reproducible micromamba environment
* Paper 1 fully runnable

