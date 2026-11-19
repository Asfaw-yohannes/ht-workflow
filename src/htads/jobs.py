"""Write VASP/SLURM job folders and submit them."""

import pprint
import subprocess
from pathlib import Path
from string import Template

DONE_FILE = "relaxed_pbe.e"      # written by the job script on success
SUBMITTED_FILE = ".submitted"    # written here after a successful sbatch

TEMPLATE_DIR = Path(__file__).resolve().parents[2] / "templates"


def load_template(name, template_dir=None):
    path = Path(template_dir or TEMPLATE_DIR) / name
    return Template(path.read_text())


def job_status(job_dir):
    """'done', 'submitted' or 'new'."""
    job_dir = Path(job_dir)
    if (job_dir / DONE_FILE).exists():
        return "done"
    if (job_dir / SUBMITTED_FILE).exists():
        return "submitted"
    return "new"


def render_slurm(template, slurm_cfg, job_name):
    """Fill the SLURM template from the `slurm` section of config.yaml."""
    account = slurm_cfg.get("account")
    modules = slurm_cfg.get("modules") or []
    values = dict(
        job_name=job_name,
        nodes=slurm_cfg["nodes"],
        ntasks_per_node=slurm_cfg["ntasks_per_node"],
        mem=slurm_cfg["mem"],
        time=slurm_cfg["time"],
        partition=slurm_cfg["partition"],
        account_line=f"#SBATCH --account={account}" if account else "",
        module_lines="\n".join(f"module load {m}" for m in modules),
        vasp_command=slurm_cfg["vasp_command"],
        vasp_pp_path=slurm_cfg["vasp_pp_path"],
        python=slurm_cfg.get("python", "python3"),
    )
    return template.substitute(values)   # substitute (not safe_) -> fails loudly


def render_vasp_script(template, structure_file, kpoints, vasp_settings, magmoms):
    return template.substitute(
        structure_file=structure_file,
        # Python literals (not JSON) so True/False stay valid Python
        kpoints=repr([int(k) for k in kpoints]),
        vasp_settings=pprint.pformat(dict(vasp_settings), indent=4, sort_dicts=False),
        magmoms=repr(dict(magmoms)),
    )


def write_job(job_dir, atoms, job_name, kpoints, cfg, slurm_tmpl, vasp_tmpl,
              vasp_settings=None):
    """
    Create one job folder: <job_name>_initial.traj, relax.py, job.slurm.
    `vasp_settings` defaults to the `vasp` section of the config (slabs).
    """
    if vasp_settings is None:
        vasp_settings = cfg["vasp"]
    job_dir = Path(job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)
    structure_file = f"{job_name}_initial.traj"
    atoms.write(job_dir / structure_file)
    (job_dir / "relax.py").write_text(
        render_vasp_script(vasp_tmpl, structure_file, kpoints,
                           vasp_settings, cfg.get("magmoms", {})))
    (job_dir / "job.slurm").write_text(render_slurm(slurm_tmpl, cfg["slurm"], job_name))


def submit(job_dir):
    """Run sbatch inside job_dir; record the job id. Raises on failure."""
    job_dir = Path(job_dir)
    res = subprocess.run(["sbatch", "--parsable", "job.slurm"], cwd=job_dir,
                         capture_output=True, text=True, check=True)
    job_id = res.stdout.strip().split(";")[0]
    (job_dir / SUBMITTED_FILE).write_text(job_id + "\n")
    return job_id
