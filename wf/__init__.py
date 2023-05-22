"""
"""

import subprocess

from latch import large_task, workflow
from latch.types import (
    LatchAuthor,
    LatchDir,
    LatchMetadata,
    LatchParameter,
)

@large_task
def rshinyA_task(
    archrproject: LatchDir,
    out_dir: str
) -> LatchDir:
    
    _r_cmd = [
        "Rscript",
        "/root/wf/task.R",
        archrproject.local_path,
    ]
    
    subprocess.run(["mkdir", "outs"])
    subprocess.run(_r_cmd)
    subprocess.run(["mv", "outs", f"{out_dir}"])

    return LatchDir(f"/root/{out_dir}", f"latch:///rshinyA_outs/{out_dir}")

metadata = LatchMetadata(
    display_name="rshiny_A",
    author=LatchAuthor(
        name="AtlasXomics Inc",
        email="jamesm@atlasxomics.com",
        github="https://github.com/atlasxomics",
    ),
    repository="https://github.com/atlasxomics/rshiny_a",
    license="MIT",
    parameters={
        "archrproject": LatchParameter(
            display_name="ArchRProject",
            description='path to ArchRProject directory in Latch, containing \
                a Save-ArchR-Project.rds file; project must have 2 conditions \
                and > 1 sample',
            batch_table_column=True, 
        ),
        "out_dir": LatchParameter(
            display_name="output directory",
            description="output directory name in rshinyA_outs/",
            batch_table_column=True,
        ),
    },
)

@workflow(metadata)
def latch_workflow(
    archrproject: LatchDir,
    out_dir: str
) -> LatchDir:
    
    return rshinyA_task(
        archrproject=archrproject,
        out_dir=out_dir
    )

if __name__ == '__main__':
    rshinyA_task(
        archrproject=LatchDir("latch:///archr_outs/craft-test2/craft-test2_25000/craft-test2_25000_ArchRProject"),
        out_dir="dev"
    )