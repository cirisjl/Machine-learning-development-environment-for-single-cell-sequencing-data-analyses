import os
import sys
import glob
import csv
import re
from pathlib import Path
from collections import defaultdict


def create_samplesheet(input_dir, output=None, prefix=None, fullpath=False, sanitize=False):
    input_dir = Path(input_dir)

    if not input_dir.is_dir():
        raise ValueError(f"Cannot read input directory: {input_dir}")

    samples = defaultdict(lambda: defaultdict(dict))

    for filename in os.listdir(input_dir):
        lower = filename.lower()

        if not lower.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz")):
            continue

        normalized = filename.replace(" ", "_")

        match = re.search(r"^(.*)_[Rr]?(1|2)(?:_|\.|$)", normalized)
        if not match:
            print(
                f"WARNING: Could not determine sample/read information from: {filename}",
                file=sys.stderr,
            )
            continue

        sample_id = match.group(1)
        read_num = match.group(2)

        lane_match = re.search(r"_L(\d{3})", sample_id, re.IGNORECASE)

        if lane_match:
            lane = lane_match.group(1)
            sample_id = re.sub(
                r"_L\d{3}",
                "",
                sample_id,
                count=1,
                flags=re.IGNORECASE,
            )
        else:
            lane = "001"

        samples[sample_id][lane][read_num] = filename

    lines = []

    for sample in sorted(samples):
        sample_name = sample.replace("//", "/", 1)

        if sanitize:
            sample_name = sample_name.replace(" ", "_")

        line = [sample_name]

        for lane in sorted(samples[sample], key=int):
            for read_num in sorted(samples[sample][lane], key=int):
                filename = samples[sample][lane][read_num]

                if fullpath:
                    filepath = str((input_dir / filename).resolve())
                elif prefix:
                    filepath = f"{prefix.rstrip('/')}/{filename}"
                else:
                    filepath = f"{str(input_dir).rstrip('/')}/{filename}"

                line.append(filepath)

        lines.append(",".join(line))

    result = "\n".join(lines)

    if result:
        result += "\n"

    if output:
        if not os.path.exists(output):
            os.makedirs(output, exist_ok=True)
        with open(output, "w", encoding="utf-8") as f:
            f.write(result)
    else:
        if not os.path.exists(output):
            os.makedirs(output, exist_ok=True)
        output = os.path.join(os.path.dirname(input_dir), 'samplesheet.csv')
        with open(output, "w", encoding="utf-8") as f:
            f.write(result)

    return output