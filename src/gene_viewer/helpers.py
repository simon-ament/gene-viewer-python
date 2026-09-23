import hashlib


def _hash_file(file_path: str, chunk_size: int = 8192) -> str:
    """Generates the SHA256 hash of a file."""

    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            sha256.update(chunk)
    return sha256.hexdigest()


def _parse_GTF_line(line: str) -> dict:
    """Parses a single line of a GTF file and returns a dictionary of its fields."""
    fields = line.strip().split("\t")
    if len(fields) != 9:
        print(f"Parsing GTF line: {line.strip()}")
        raise ValueError("Invalid GTF line: must have 9 fields.")

    return {
        "seqname": fields[0],
        "source": fields[1],
        "feature": fields[2],
        "start": int(fields[3]),
        "end": int(fields[4]),
        "score": fields[5] if fields[5] != "." else None,
        "strand": fields[6] if fields[6] != "." else None,
        "frame": fields[7] if fields[7] != "." else None,
        "attributes": fields[8] if fields[8] != "." else None,
    }


def _get_GTF_attribute(attributes: str, key: str):
    # fast search for the key in the attributes string (not using regex for performance reasons)
    key_pattern = f'{key} "'
    start_index = attributes.find(key_pattern)
    if start_index == -1:
        return None
    # Find the end of the value (the closing quote)
    end_index = attributes.find('"', start_index + len(key_pattern))
    if end_index == -1:
        return None
    return attributes[start_index + len(key_pattern) : end_index]
