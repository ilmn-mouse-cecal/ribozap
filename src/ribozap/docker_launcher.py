import subprocess
import logging
from pathlib import Path
from os import path, environ

logger = logging.getLogger(__name__)


def run_docker(image: str, mount_file: Path, out_dir: Path, container_cmd: str, cpus=4, memory=16):

    docker_options = []

    if cpus:
        docker_options.append(f"--cpus={cpus}")
    if memory:
        docker_options.append(f"--memory={memory}g")

    # Read mount lines from file
    try:
        with mount_file.open() as f:
            mount_lines = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        logger.error(f"Mount file not found: {mount_file}")
        raise

    # Make sure /work directory exists
    project_root = Path(__file__).resolve().parents[2]
    nextflow_files = project_root / "nextflow_files"
    nextflow_logs = project_root / "logs"

    # Append -v mounts for /app
    mount_lines.append(f'-v "{out_dir.resolve()}:/app/{out_dir.name}"')
    mount_lines.append(f'-v "{nextflow_files.resolve()}:/app/nextflow_files"')
    mount_lines.append(f'-v "{nextflow_logs.resolve()}:/app/logs"')
    

    # Add rRNA index
    rrna_index = Path(environ["HOME"]) / ".rRNA_database_index"
    mount_lines.append(f'-v "{rrna_index}:/app/idx"')

    # Final Docker command
    docker_cmd = (
        f'docker run --rm {" ".join(docker_options)} {" ".join(mount_lines)} '
        f'{image} bash -c "cd /app && {container_cmd}"'
    )

    logger.info("🐳 Running Docker command:")
    logger.info(docker_cmd)

    try:
        subprocess.run(docker_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Docker command failed: {e}")
        raise
