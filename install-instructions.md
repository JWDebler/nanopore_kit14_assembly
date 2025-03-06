# Installation Guide: Conda, Mamba, Dorado, and Basecalling Environment

This guide provides step-by-step instructions to install **Conda, Mamba, Dorado**, and set up a **basecalling environment**.

---

## **:one: nstall Conda (Miniconda)**
Miniconda is a lightweight version of Anaconda that includes only Conda.

### **For Linux/macOS:**
1. **Download Miniconda**:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   ```
2. **Run the installer**:
   ```bash
   bash Miniconda3-latest-Linux-x86_64.sh
   ```
   - Press **Enter** to continue.
   - Accept the license (`yes`).
   - Choose an **installation path** (default: `$HOME/miniconda3`).
   - Allow the installer to initialize Conda in your shell.

3. **Restart your shell**:
   ```bash
   source ~/.bashrc
   ```
4. **Verify installation**:
   ```bash
   conda --version
   ```

---

## **:two: Install Mamba (Faster Conda)**
Mamba is a fast, drop-in replacement for Conda.

1. **Install Mamba**:
   ```bash
   conda install -n base -c conda-forge mamba
   ```
2. **Initialize Mamba** in your shell:
   ```bash
   eval "$(mamba shell hook --shell bash)"
   mamba shell init --shell bash --root-prefix=~/.local/share/mamba
   source ~/.bashrc
   ```
3. **Verify Mamba installation**:
   ```bash
   mamba --version
   ```

---

## ** :three: Install Dorado**
Dorado is an Oxford Nanopore basecalling tool.

1. **Download Dorado**:
   ```bash
   wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz
   ```
2. **Move Dorado to a system-wide location**:
   ```bash
   sudo mv dorado-0.9.1-linux-x64 /opt/dorado
   ```
3. **Add Dorado to your PATH**:
   ```bash
   echo 'export PATH=/opt/dorado/bin/:$PATH' >> ~/.bashrc
   source ~/.bashrc
   ```
4. **Verify Dorado installation**:
   ```bash
   dorado --version
   ```

---

## ** :four: reate and Activate the Conda Environment**
We will now set up a **Conda environment** using Mamba.

### **Create the environment**
1. Ensure you have an `environment.yml` file with the following contents:

   ```yaml
   name: basecalling
   channels:
     - bioconda
     - conda-forge
     - jannessp
     - defaults
   dependencies:
     - bioconda::samtools=1.21
     - jannessp::pod5=0.3.23
     - bioconda::chopper=0.9.1
     - conda-forge::pigz=2.8
     - pip
   variables:
     PATH: "/opt/dorado/bin:$PATH"
   ```

2. **Create the environment using Mamba**:
   ```bash
   mamba env create -f environment.yml
   ```

3. **Activate the environment**:
   ```bash
   mamba activate basecalling
   ```

---

## **:five: Verify Installation**
After installing all components, perform the following checks:

1. **Check installed Conda packages**:
   ```bash
   conda list
   ```
   This should display all installed packages, including `samtools`, `chopper`, `pod5`, and `pigz`.

2. **Check if Dorado is accessible**:
   ```bash
   dorado --version
   ```

3. **Check if Mamba can manage the environment**:
   ```bash
   mamba list
   ```

---

## **:six: Summary of Installed Tools**
| Tool   | Installation Method |
|--------|---------------------|
| **Conda (Miniconda)** | Installed from Miniconda script |
| **Mamba** | Installed via Conda (`conda install -c conda-forge mamba`) |
| **Dorado** | Downloaded and manually added to PATH |
| **Basecalling Environment** | Created with Mamba (`mamba env create -f environment.yml`) |

---

## **Next Steps**
Now that your environment is set up, you can start **basecalling** with Dorado or other tools. If you encounter any issues, ensure that:
- Mamba is correctly initialized (`eval "$(mamba shell hook --shell bash)"`).
- Dorado is correctly linked (`echo $PATH` should include `/opt/dorado/bin`).
- The correct Conda environment is activated (`mamba activate basecalling`).
