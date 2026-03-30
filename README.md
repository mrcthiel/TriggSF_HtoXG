g TriggSF HtoXG

Trigger scale factor for **Mu17-Pho30** using NanoAOD and coffea.

---

## 1. Clone repository

```bash
git clone https://github.com/mrcthiel/TriggSF_HtoXG.git
cd TriggSF_HtoXG
```

---

## 2. Setup environment

This project uses **micromamba**.

### Create environment
```bash
micromamba create -f environment.yml
```

### Activate
```bash
micromamba activate hztest
```

---

## 3. Run analysis

### Full workflow
```bash
python run_analysis_trigg.py
```

