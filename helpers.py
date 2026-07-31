import hashlib


def _hash_file(file_path: str, chunk_size: int = 8192) -> str:
    """Generates the SHA256 hash of a file."""

    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        while chunk := f.read(chunk_size):
            sha256.update(chunk)
    return sha256.hexdigest()
