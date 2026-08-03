#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import json
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

NOTEBOOK_APPS: Mapping[str, tuple[str, ...]] = {
    "01_prerequisites.ipynb": (),
    "02_initialization.ipynb": ("barebone_template",),
    "03_transformation.ipynb": ("backward_trans", "forward_trans"),
    "04_operation.ipynb": ("laplacian", "inverse_laplacian"),
    "05_time_integration.ipynb": (
        "time_integration_first",
        "time_integration_second",
    ),
    "06_vector_field.ipynb": ("vecfld_reconstruction", "tp_project"),
}

DOCUMENT_APPS: Mapping[str, tuple[str, ...]] = {
    "docs/release_notes_v1.1.x.md": (),
    "docs/release_notes_v1.0.x.md": (),
    "docs/tutorial/prerequisites.md": (),
    "docs/sample_programs/1d_radial_wave_propagation.md": (
        "wave_propagation_1d",
    ),
    "docs/sample_programs/2d_scalar_diffusion.md": ("scalar_diffusion_2d",),
    "docs/sample_programs/2d_scalar_transport.md": ("scalar_transport_2d",),
    "docs/sample_programs/3d_vortical_flow.md": ("vortical_flow_3d",),
    "docs/tutorial/initialization.md": ("barebone_template",),
    "docs/tutorial/transformation.md": ("backward_trans", "forward_trans"),
    "docs/tutorial/operation.md": ("laplacian", "inverse_laplacian"),
    "docs/tutorial/time_integration.md": (
        "time_integration_first",
        "time_integration_second",
    ),
    "docs/tutorial/vector_field.md": ("vecfld_reconstruction", "tp_project"),
}

ALL_APPS = tuple(dict.fromkeys(app for apps in DOCUMENT_APPS.values() for app in apps))

# Applications that must be launched on exactly one process, whatever --ranks says.
# wave_propagation_1d solves a purely radial problem, so the rank that owns the radial
# direction has to own all of it; the program stops deliberately when given more than one
# process. Running the whole gate at --ranks 1 to accommodate it would give up parallel
# coverage of the other twelve, so the rank count is overridden per application instead.
SERIAL_ONLY_APPS = frozenset({"wave_propagation_1d"})


class ValidationError(RuntimeError):
    pass


def _cell_source(cell: Mapping[str, object]) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(str(line) for line in source)
    return str(source)


def validate_notebooks(root: Path) -> None:
    tutorial_dir = root / "tutorials"
    found = {path.name for path in tutorial_dir.glob("*.ipynb")}
    expected = set(NOTEBOOK_APPS)
    missing = sorted(expected - found)
    unexpected = sorted(found - expected)
    if missing or unexpected:
        raise ValidationError(
            f"notebook inventory mismatch; missing={missing}, unexpected={unexpected}"
        )

    for name, apps in NOTEBOOK_APPS.items():
        path = tutorial_dir / name
        try:
            notebook = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValidationError(f"cannot read {path}: {exc}") from exc

        cells = notebook.get("cells")
        if not isinstance(cells, list):
            raise ValidationError(f"{path} has no notebook cell list")

        sources = [_cell_source(cell) for cell in cells if isinstance(cell, dict)]
        combined = "\n".join(sources)
        stale = [token for token in ("mpirun.openmpi", "nproc", "--oversubscribe") if token in combined]
        if stale:
            raise ValidationError(f"{path} contains stale platform command(s): {stale}")

        for app in apps:
            if f"make {app}" not in combined or f"build/bin/{app}" not in combined:
                raise ValidationError(
                    f"{path} does not contain both build and run commands for {app}"
                )
        if name == "01_prerequisites.ipynb" and "make mods" not in combined:
            raise ValidationError(f"{path} does not contain the module build command")

        for cell in cells:
            if not isinstance(cell, dict) or cell.get("cell_type") != "code":
                continue
            source = _cell_source(cell)
            if source.lstrip().startswith("%%bash"):
                continue
            try:
                ast.parse(source, filename=str(path))
            except SyntaxError as exc:
                raise ValidationError(f"Python syntax error in {path}: {exc}") from exc


def validate_documents(root: Path) -> None:
    for relative, apps in DOCUMENT_APPS.items():
        path = root / relative
        if not path.is_file():
            raise ValidationError(f"missing documentation page: {relative}")
        text = path.read_text(encoding="utf-8")
        stale = [token for token in ("mpirun.openmpi", "nproc", "--oversubscribe") if token in text]
        if stale:
            raise ValidationError(f"{relative} contains stale platform command(s): {stale}")
        for app in apps:
            if f"make {app}" not in text or f"build/bin/{app}" not in text:
                raise ValidationError(
                    f"{relative} does not contain both build and run commands for {app}"
                )


def _run(
    command: Sequence[str],
    *,
    cwd: Path,
    timeout: int,
    label: str,
    env: Mapping[str, str] | None = None,
) -> str:
    try:
        completed = subprocess.run(
            list(command),
            cwd=cwd,
            env=dict(env) if env is not None else None,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ValidationError(f"{label} failed to start or timed out: {exc}") from exc
    if completed.returncode:
        output = completed.stdout.strip().splitlines()
        tail = "\n".join(output[-30:])
        raise ValidationError(f"{label} exited {completed.returncode}:\n{tail}")
    return completed.stdout


def _input_text(
    *,
    nr: int,
    np: int,
    nz: int,
    nrchop: int,
    npchop: int,
    nzchop: int,
    ell: float,
    dt: float,
    totaltime: float,
    totaln: int,
    visc: float,
    hyperpow: int = 0,
    hypervisc: float = 0.0,
    fields: bool = False,
    field_interval: int = 1,
    logs: bool = False,
    log_interval: int = 1,
    stability: bool = True,
    svv_cutoff: float = 0.75,
    svv_target: float = 2.0e-2,
) -> str:
    return f"""!!! COMPUTATIONAL DOMAIN INFO !!!
# ---------- NR ----------- NP ----------- NZ ----------------------------------
             {nr}             {np}              {nz}
# ------ NRCHOP ------- NPCHOP ------- NZCHOP ----------------------------------
             {nrchop}             {npchop}              {nzchop}
# --------- ELL --------- ZLEN ------ ZLENxPI ---(IF ZLENxPI==T, ZLEN=ZLEN*PI)--
           {ell:.6E}           1.D0              F
#
!!! TIME STEPPING INFO !!!
# ---------- DT ----------- TI --------- TOTT ----------- NI --------- TOTN ----
          {dt:.6E}           0.D0          {totaltime:.6E}              0             {totaln}
#
!!! FIELD PROPERTY INFO !!!
# -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
          {visc:.6E}              {hyperpow}           {hypervisc:.6E}
#
!!! FIELD SAVE INFO !!!
# --------------------- FLDDIR -------------------------------------------------
                  ./output/fld
# ---- ISFLDSAV -- FLDSAVINTVL ---(IF ISFLDSAV!=T, FIELDS ARE NOT SAVED)--------
              {'T' if fields else 'F'}              {field_interval}
#
!!! DATA LOGGING INFO !!!
# --------------------- LOGDIR -------------------------------------------------
                  ./output/dat
# ---- ISLOGSAV -- LOGSAVINTVL ---(IF ISLOGSAV!=T, LOGS ARE NOT GENERATED)------
              {'T' if logs else 'F'}              {log_interval}
!!! STABILITY CONTROL INFO (OPTIONAL) !!!
# ------- IS_SVV ------- CUTOFF ------- TARGET ------ STRENGTH ------- RELAX ---
              {'T' if stability else 'F'}          {svv_cutoff:.6E}          {svv_target:.6E}          0.12          0.25
/* ------------------------------ END OF INPUT ------------------------------ */
"""


def _prepare_run_directory() -> Path:
    run_dir = Path(tempfile.mkdtemp(prefix="mlegs-tutorial-"))
    (run_dir / "output" / "fld").mkdir(parents=True)
    (run_dir / "output" / "dat").mkdir(parents=True)

    configs = {
        "input.params": _input_text(
            nr=32,
            np=16,
            nz=8,
            nrchop=32,
            npchop=9,
            nzchop=5,
            ell=4.0,
            dt=1.0e-2,
            totaltime=2.0e-2,
            totaln=2,
            visc=1.0e-4,
            hyperpow=8,
            hypervisc=5.0e-7,
            fields=False,
            logs=False,
        ),
        "input_1d.params": _input_text(
            nr=32,
            np=1,
            nz=1,
            nrchop=32,
            npchop=1,
            nzchop=1,
            ell=1.0,
            dt=1.0e-2,
            totaltime=2.0e-2,
            totaln=2,
            visc=0.0,
            fields=False,
            logs=False,
        ),
        "input_2d.params": _input_text(
            nr=32,
            np=48,
            nz=1,
            nrchop=32,
            npchop=25,
            nzchop=1,
            ell=1.0,
            dt=1.0e-2,
            totaltime=2.0e-2,
            totaln=2,
            visc=5.0e-3,
            fields=False,
            logs=False,
        ),
        "input_tutorial.params": _input_text(
            nr=32,
            np=48,
            nz=1,
            nrchop=32,
            npchop=25,
            nzchop=1,
            ell=1.0,
            dt=2.0e-1,
            totaltime=4.0e-1,
            totaln=2,
            visc=1.0e-1,
            fields=True,
            field_interval=1,
            logs=False,
            stability=False,
        ),
        "input_svv_off.params": _input_text(
            nr=32,
            np=16,
            nz=8,
            nrchop=32,
            npchop=9,
            nzchop=5,
            ell=4.0,
            dt=5.0e-3,
            totaltime=1.0e-1,
            totaln=20,
            visc=1.0e-10,
            hyperpow=8,
            hypervisc=1.0e-10,
            fields=True,
            field_interval=20,
            logs=False,
            stability=False,
            svv_cutoff=0.5,
            svv_target=1.0e-12,
        ),
        "input_svv_on.params": _input_text(
            nr=32,
            np=16,
            nz=8,
            nrchop=32,
            npchop=9,
            nzchop=5,
            ell=4.0,
            dt=5.0e-3,
            totaltime=1.0e-1,
            totaln=20,
            visc=1.0e-10,
            hyperpow=8,
            hypervisc=1.0e-10,
            fields=True,
            field_interval=20,
            logs=False,
            stability=True,
            svv_cutoff=0.5,
            svv_target=1.0e-12,
        ),
    }
    for name, text in configs.items():
        path = run_dir / name
        path.write_text(text, encoding="utf-8")
    return run_dir


def _mpi_command(mpirun: str, ranks: int) -> list[str]:
    command = [mpirun]
    version = subprocess.run(
        [mpirun, "--version"],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    ).stdout.lower()
    if "open mpi" in version or "openmpi" in version:
        command.append("--oversubscribe")
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            command.append("--allow-run-as-root")
    command.extend(["-n", str(ranks)])
    return command


def run_executables(root: Path, *, ranks: int, timeout: int, skip_build: bool) -> None:
    mpirun = os.environ.get("MPIEXEC") or shutil.which("mpirun")
    if not mpirun:
        raise ValidationError("mpirun is required for executable tutorial validation")

    if not skip_build:
        _run(
            ["make", *ALL_APPS],
            cwd=root,
            timeout=timeout,
            label="tutorial executable build",
        )

    binaries = {app: root / "build" / "bin" / app for app in ALL_APPS}
    missing = [app for app, path in binaries.items() if not path.is_file()]
    if missing:
        raise ValidationError(
            "missing tutorial executable(s): "
            + ", ".join(missing)
            + "; run external/CMake_build.sh and make first"
        )

    run_dir = _prepare_run_directory()
    try:
        mpi = _mpi_command(mpirun, ranks)

        serial = _mpi_command(mpirun, 1) if ranks != 1 else mpi

        def run_case(app: str, config: str) -> None:
            if not (run_dir / config).is_file():
                raise ValidationError(f"missing input configuration: {config}")
            launcher = serial if app in SERIAL_ONLY_APPS else mpi
            _run(
                [*launcher, str(binaries[app])],
                cwd=run_dir,
                timeout=timeout,
                label=f"{app} tutorial",
                env=os.environ,
            )

        run_case("barebone_template", "input.params")
        run_case("backward_trans", "input_2d.params")
        run_case("forward_trans", "input_2d.params")
        run_case("laplacian", "input_2d.params")
        run_case("inverse_laplacian", "input_2d.params")
        run_case("time_integration_first", "input_tutorial.params")
        run_case("time_integration_second", "input_tutorial.params")

        vector_input = _input_text(
            nr=200,
            np=128,
            nz=8,
            nrchop=200,
            npchop=65,
            nzchop=5,
            ell=5.0,
            dt=1.0e-1,
            totaltime=1.0e-1,
            totaln=1,
            visc=1.0e-3,
            fields=True,
            field_interval=1,
            logs=False,
            stability=False,
        )
        (run_dir / "input_tutorial.params").write_text(vector_input, encoding="utf-8")
        run_case("vecfld_reconstruction", "input_tutorial.params")
        run_case("tp_project", "input_tutorial.params")

        run_case("wave_propagation_1d", "input_1d.params")
        run_case("scalar_diffusion_2d", "input_2d.params")
        run_case("scalar_transport_2d", "input_2d.params")
        run_case("vortical_flow_3d", "input.params")

        svv_states = []
        svv_field = run_dir / "output" / "fld" / "vormagfld000020_PPP"
        for config in ("input_svv_off.params", "input_svv_on.params"):
            shutil.copyfile(run_dir / config, run_dir / "input.params")
            run_case("vortical_flow_3d", "input.params")
            if not svv_field.is_file():
                raise ValidationError(f"missing SVV validation field: {svv_field.name}")
            svv_states.append(svv_field.read_bytes())
        if svv_states[0] == svv_states[1]:
            raise ValidationError("SVV enabled/disabled validation fields are identical")

        rejected_inputs = {
            "zero enabled hyperviscosity": _input_text(
                nr=32, np=16, nz=8, nrchop=32, npchop=9, nzchop=5,
                ell=4.0, dt=1.0e-2, totaltime=1.0e-2, totaln=1,
                visc=1.0e-4, hyperpow=8, hypervisc=0.0,
            ),
            "unsupported hyperviscosity power": _input_text(
                nr=32, np=16, nz=8, nrchop=32, npchop=9, nzchop=5,
                ell=4.0, dt=1.0e-2, totaltime=1.0e-2, totaln=1,
                visc=1.0e-4, hyperpow=10, hypervisc=1.0e-7,
            ),
            "invalid SVV cutoff": _input_text(
                nr=32, np=16, nz=8, nrchop=32, npchop=9, nzchop=5,
                ell=4.0, dt=1.0e-2, totaltime=1.0e-2, totaln=1,
                visc=1.0e-4, hyperpow=8, hypervisc=1.0e-7,
                svv_cutoff=2.0,
            ),
        }
        for label, text in rejected_inputs.items():
            (run_dir / "input.params").write_text(text, encoding="utf-8")
            try:
                completed = subprocess.run(
                    [*mpi, str(binaries["vortical_flow_3d"])],
                    cwd=run_dir,
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=timeout,
                )
            except (OSError, subprocess.TimeoutExpired) as exc:
                raise ValidationError(f"{label} rejection check failed: {exc}") from exc
            if completed.returncode == 0:
                raise ValidationError(f"{label} was accepted")
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def execute_notebooks(root: Path, *, timeout: int) -> None:
    try:
        import nbformat
        from nbclient import NotebookClient
    except ImportError as exc:
        raise ValidationError(
            "--execute-notebooks requires nbformat and nbclient; install them "
            "alongside the notebook plotting dependencies"
        ) from exc

    for name in NOTEBOOK_APPS:
        path = root / "tutorials" / name
        notebook = nbformat.read(path, as_version=4)
        client = NotebookClient(
            notebook,
            timeout=timeout,
            kernel_name="python3",
            resources={"metadata": {"path": str(root / "tutorials")}},
        )
        try:
            client.execute()
        except Exception as exc:
            raise ValidationError(f"notebook execution failed for {name}: {exc}") from exc


def parse_args(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate MLegS tutorials")
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="reuse existing build/bin executables instead of invoking make",
    )
    parser.add_argument(
        "--execute-notebooks",
        action="store_true",
        help="also execute notebook cells (requires nbformat and nbclient)",
    )
    parser.add_argument(
        "--ranks",
        type=int,
        default=2,
        help="MPI ranks for executable smoke tests (default: 2)",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="timeout in seconds for each build/run operation (default: 600)",
    )
    return parser.parse_args(list(argv))


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.ranks < 1 or args.timeout < 1:
        print("--ranks and --timeout must be positive", file=sys.stderr)
        return 2
    root = Path(__file__).resolve().parents[1]
    try:
        validate_notebooks(root)
        validate_documents(root)
        if args.execute_notebooks:
            execute_notebooks(root, timeout=args.timeout)
        run_executables(
            root,
            ranks=args.ranks,
            timeout=args.timeout,
            skip_build=args.skip_build,
        )
    except ValidationError as exc:
        print(f"tutorial validation failed: {exc}", file=sys.stderr)
        return 1
    print(
        f"tutorial validation passed: {len(NOTEBOOK_APPS)} notebooks, "
        f"{len(DOCUMENT_APPS)} documentation pages, "
        f"and {len(ALL_APPS)} executable tutorial paths"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
