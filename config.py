import os


DEFAULT_CONFIG_FILE = os.path.join(
    "example",
    "amatrice",
    "input_file.txt",
)


def _parse_value(value):
    value = value.strip()
    if value == "None":
        return None
    return value


def _optional_float(value):
    if value is None:
        return None
    return float(value)


def parse_int_list(value, default=None):
    """
    Parse comma-separated integers and inclusive ranges.

    Examples
    --------
    "1,3,5-7" -> [1, 3, 5, 6, 7]

    Duplicate values are removed while preserving their first occurrence.
    """
    if value is None:
        return default

    parsed = []

    for chunk in str(value).split(","):
        chunk = chunk.strip()

        if not chunk:
            continue

        if "-" in chunk:
            start_text, end_text = chunk.split("-", 1)
            start = int(start_text)
            end = int(end_text)

            if end < start:
                raise ValueError(
                    f"Invalid descending range: {chunk}"
                )

            parsed.extend(range(start, end + 1))

        else:
            parsed.append(int(chunk))

    if not parsed:
        raise ValueError(
            f"No integer values found in: {value!r}"
        )

    # Remove duplicates while preserving order
    return list(dict.fromkeys(parsed))


def load_config(file=DEFAULT_CONFIG_FILE):
    """
    Read and validate a RESTORE input file.
    """
    config_file = os.path.abspath(file)
    config_dir = os.path.dirname(config_file)

    allowed_keys = {
        "name_catalog",
        "delimiter",
        "sizes",
        "steps",
        "st_dev_multiplier",
        "Sigma",
        "sbin",
        "fault_length_multiplier",
        "t_end_quiet",
        "b-value",
        "alpha",
        "mc",
        "depth_distribution",
        "p0",
        "seeds",
        "output_folder",
    }

    values = {}

    with open(config_file, "r") as input_file:
        for raw_line in input_file:
            line = raw_line.strip()

            if (
                not line
                or line.startswith("*")
                or line.startswith("-")
                or ":" not in line
            ):
                continue

            key, value = line.split(":", 1)
            key = key.strip()

            if key not in allowed_keys:
                raise ValueError(
                    f"Unknown configuration key: {key}"
                )

            values[key] = _parse_value(value)

    required = {
        "name_catalog",
        "delimiter",
        "sizes",
        "steps",
        "st_dev_multiplier",
        "Sigma",
        "sbin",
        "fault_length_multiplier",
        "t_end_quiet",
        "alpha",
        "depth_distribution",
        "p0",
    }

    missing = sorted(required - values.keys())

    if missing:
        raise ValueError(
            "Missing required configuration keys: "
            + ", ".join(missing)
        )

    # Resolve catalog path relative to the input file
    name_catalog = str(values["name_catalog"])

    if not os.path.isabs(name_catalog):
        name_catalog = os.path.abspath(
            os.path.join(config_dir, name_catalog)
        )

    # Resolve output path relative to the input file
    output_folder = str(
        values.get(
            "output_folder",
            os.path.join(
                "..",
                "..",
                "outputs",
                "restore_run",
            ),
        )
    )

    if not os.path.isabs(output_folder):
        output_folder = os.path.abspath(
            os.path.join(config_dir, output_folder)
        )

    # Parameter values to explore
    sizes = parse_int_list(values["sizes"])
    steps = parse_int_list(values["steps"])
    seeds = parse_int_list(
        values.get("seeds"),
        default=[1],
    )

    if any(size < 50 for size in sizes):
        raise ValueError(
            "Every window size must be at least 50 events."
        )

    if any(step <= 0 for step in steps):
        raise ValueError(
            "Every window step must be positive."
        )

    return {
        "config_file": config_file,
        "name_catalog": name_catalog,
        "delimiter": str(values["delimiter"]),
        "sizes": sizes,
        "steps": steps,
        "seeds": seeds,
        "st_dev_multiplier": float(
            values["st_dev_multiplier"]
        ),
        "Sigma": float(values["Sigma"]),
        "sbin": float(values["sbin"]),
        "fault_length_multiplier": int(
            values["fault_length_multiplier"]
        ),
        "t_end_quiet": str(values["t_end_quiet"]),
        "b-value": _optional_float(values.get("b-value")),
        "alpha": float(values["alpha"]),
        "mc": _optional_float(values.get("mc")),
        "depth_distribution": str(
            values["depth_distribution"]
        ),
        "p0": [
            float(value)
            for value in str(values["p0"]).split(",")
        ],
        "output_folder": output_folder,
    }
