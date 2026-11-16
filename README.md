# vizbench-r

R (and helper Python) code for the **vizbench** benchmark of single-cell normalization → integration → visualization pipelines.

```
vizbench-r/
├── viz_benchmark.yaml   # omnibenchmark plan
├── run_vizbench.R       # single entry point (--what / --flavour)
├── utils/               # <what>_utils.R|py, one per stage
├── envs/                # Apptainer definition (+ placeholder conda yml)
└── explorer/            # index.html results explorer + pipeline_data.json
```

## 1. Build the environment

This benchmark uses the **Apptainer** backend and expects the image at `~/envs/r-451-bioc-321-apptainer.sif`

```bash
git clone https://github.com/zhiqianZ/vizbench-r.git
cd vizbench-r

mkdir -p ~/envs
cp envs/r-451-bioc-321-dummy.yml ~/envs/
apptainer build ~/envs/r-451-bioc-321-apptainer.sif envs/r-451-bioc-321-apptainer.def
```

## 2. Install omnibenchmark

Follow the official installation guide: **<https://docs.omnibenchmark.org/latest>** (step-by-step [tutorial](https://omnibenchmark.org/tutorial/)). In short:

```bash
pip install omnibenchmark      # Python env; miniforge recommended
ob --version
```

## 3. Run the benchmark

From the repo root:

```bash
nohup ob run benchmark -b viz_benchmark.yaml --cores 10 --task-timeout 86400 &
```

Results are written to `out/`, one folder per stage/module/parameter set.

## 4. Explore the results

The aggregated results are stored in `explorer/pipeline_data.json` and can be explored from https://zhiqianz.github.io/vizbench-r/explorer/.
